import glob
import json
import numpy as np

from collections import defaultdict
from scipy.sparse import coo_matrix

from melon.utils import *
from .utils import *


class AntibioticResistanceGeneProfiler:
    def __init__(self, file, db, output, threads=os.cpu_count()):
        '''
        Initialization.
        '''
        self.file = file
        self.db = db
        self.output = output
        self.threads = threads
        self.outfile = get_filename(self.file, self.output)

        ## temporary variables for diamond's hits or minimap2's alignments
        self.hits, self.alignments = [], []

        ## taxonomy assignments
        self.assignments = {}

        ## genome copies
        self.lineage2genome = {}
        with open(f'{self.outfile}.tsv') as f:
            next(f)
            for line in f:
                ls = line.rstrip().split('\t')
                self.lineage2genome[';'.join(ls[:7])] = float(ls[7])

        ## neglect non-prokaryotic reads if necessary
        with open(f'{self.outfile}.json') as f:
            self.nset = {key for key, val in json.load(f).items() if val['remark'] == 'putatively non-prokaryotic'}

    def pre_overlap(self):
        '''
        Get overlaps of top 10k reads.
        '''
        subprocess.run([
            'seqkit', 'head',
            '-n', '10000',
            '-o', f'{self.outfile}.init.sequence.tmp',
            self.file
        ], check=True, stderr=subprocess.DEVNULL)

        subprocess.run([
            'minimap2',
            '-x', 'ava-ont',
            '-t', str(self.threads),
            '-o', f'{self.outfile}.init.overlap.tmp',
            f'{self.outfile}.init.sequence.tmp', f'{self.outfile}.init.sequence.tmp'
        ], check=True, stderr=subprocess.DEVNULL)

    def post_overlap(self, chunk_size=0):
        '''
        Get overlaps of ARG-containing reads.
        '''
        qseqid = {hit[0] for hit in self.hits}
        nchunks = len(qseqid) // chunk_size + 1 if chunk_size != 0 else 1

        extract_sequence(self.file, qseqid, f'{self.outfile}.sarg.sequence.tmp')
        if nchunks > 1:
            for file in glob.glob(f'{self.outfile}.sarg.sequence.part*.tmp'):
                os.remove(file)

            subprocess.run([
                'seqkit', 'split',
                '-p', str(nchunks),
                '--out-dir', self.output,
                f'{self.outfile}.sarg.sequence.tmp'
            ], check=True, stderr=subprocess.DEVNULL)

        for file in glob.glob(f'{self.outfile}.sarg.sequence.part*.tmp') if nchunks > 1 else [f'{self.outfile}.sarg.sequence.tmp']:
            subprocess.run([
                'minimap2',
                '-x', 'ava-ont',
                '-t', str(self.threads),
                '-o', file.replace('.sarg.sequence.', '.sarg.overlap.'),
                file, file
            ], check=True, stderr=subprocess.DEVNULL)

    def run_overlap(self, mode, identity=0, chunk_size=0, scale=2.5):
        '''
        Run read overlapping.
        '''
        ## initial cutoff for pre-screening
        if mode == 'pre':
            if identity == 0:
                if not self.debug: self.pre_overlap()
                overlaps = filter_overlap(f'{self.outfile}.init.overlap.tmp')

                ## default initial identity cutoff: 0.9 * (90 - 2.5 * median sequence divergence %)
                DV = np.median(np.fromiter((x[-1] for x in overlaps), dtype=float))
                self.identity = 0.9 * 100 * (0.9 - scale * DV)
                logger.info(
                    f'... median sequence divergence: {DV:.4f}'
                    f' | initial identity cutoff: 0.9 * {self.identity / 0.9:.2f}.'
                )
            else:
                self.identity = identity

        ## compute an overlap identity based on all ARG-containing reads
        if mode == 'post':
            if not self.debug: self.post_overlap(chunk_size=chunk_size)
            self.overlaps = []
            for file in glob.glob(f'{self.outfile}.sarg.overlap*.tmp'):
                self.overlaps.extend(filter_overlap(file))

            DV = np.median(np.fromiter((x[-1] for x in self.overlaps), dtype=float))
            self.DV = max(scale * DV, 0.05)

            if identity == 0:
                ## update identity cutoff if necessary
                self.identity = max(self.identity, 100 * (0.9 - scale * DV))
                nhits, self.hits = len(self.hits), [hit for hit in self.hits if hit[8] >= self.identity]
                logger.info(
                    f'... median sequence divergence of ARG-containing reads: {DV:.4f}'
                    f' | identity cutoff: {self.identity:.2f}'
                    f' | low-identity HSPs: {nhits - len(self.hits)}.'
                )
            else:
                self.identity = identity

    def run_diamond(self, max_target_seqs=25, evalue=1e-5, identity=0):
        '''
        Run diamond blastx.
        '''
        outfmt = ['qseqid', 'sseqid', 'pident', 'qlen', 'qstart', 'qend', 'slen', 'sstart', 'send', 'evalue', 'bitscore']
        subprocess.run([
            'diamond', 'blastx',
            '--db', os.path.join(self.db, 'sarg.dmnd'),
            '--query', self.file,
            '--out', f'{self.outfile}.sarg.diamond.tmp',
            '--outfmt', '6', *outfmt,
            '--evalue', str(evalue),
            '--id', str(identity),
            '--range-culling',
            '--frameshift', '15',
            '--range-cover', '25',
            '--max-hsps', '0',
            '--max-target-seqs', str(max_target_seqs),
            '--threads', str(self.threads)
        ], check=True, stderr=subprocess.DEVNULL)

    def parse_diamond(self):
        '''
        Parse diamond's output and record the hits.
        '''
        with open(f'{self.outfile}.sarg.diamond.tmp') as f:
            for line in f:
                ls = line.rstrip().split('\t')
                qseqid, sseqid = ls[0], ls[1]
                if qseqid not in self.nset:
                    qstart, qend = sort_coordinate(int(ls[4]), int(ls[5]))
                    sstart, send = sort_coordinate(int(ls[7]), int(ls[8]))
                    qlen, slen = int(ls[3]), int(ls[6])
                    scov = (send - sstart) / slen
                    bitscore = float(ls[10])

                    ## append qseqid and coordinates for back-tracing
                    if (pident := float(ls[2])) >= self.identity:
                        self.hits.append([qseqid, sseqid, qlen, qstart, qend, slen, sstart, send, pident, scov, bitscore])

        logger.info(
            f'... candidate HSPs: {len(self.hits)}'
            f' | ARG-containing reads: {len({hit[0] for hit in self.hits})}.'
        )

    def run_minimap(self, secondary_num=2147483647, secondary_ratio=0.9):
        '''
        Run minimap2 to get taxonomic profiles.
        '''
        qseqids = defaultdict(set)
        for hit in self.hits:
            qseqids[hit[1].split('|')[1]].add(hit[0])

        with open(f'{self.outfile}.sarg.minimap.tmp', 'w') as w:
            with logging_redirect_tqdm():
                bar_format = '==> {desc}{percentage:3.0f}%|{bar}|{n_fmt}/{total_fmt} [{elapsed}<{remaining}]'
                queue = tqdm(dict(sorted(qseqids.items())).items(), bar_format=bar_format, leave=False)
                for family, qseqid in queue:
                    queue.set_description(f'Processing <{family}>')
                    sequences = extract_sequence(f'{self.outfile}.sarg.sequence.tmp', qseqid)
                    subprocess.run([
                        'minimap2',
                        '-cx', 'map-ont',
                        '-f', '0',
                        '-N', str(secondary_num),
                        '-p', str(secondary_ratio),
                        '-t', str(self.threads),
                        os.path.join(self.db, f'sarg.{family}.mmi'), '-',
                    ], check=True, stdout=w, stderr=subprocess.DEVNULL, input=sequences, text=True)

    def parse_minimap(self, plasmid=False):
        '''
        Parse minimap's paf output.
        '''
        accession2lineage = {}
        with open(os.path.join(self.db, 'sarg.metadata.tsv')) as f:
            next(f)
            for line in f:
                ls = line.rstrip().split('\t')
                accession2lineage[ls[0]] = ';'.join(ls[1:])

        qcoords = defaultdict(set)
        for hit in self.hits:
            qcoords[hit[0]].add(tuple(hit[3:5]))

        alignments = []
        scores, max_scores = defaultdict(dict), {}
        with open(f'{self.outfile}.sarg.minimap.tmp') as f:
            for line in f:
                ls = line.rstrip().split('\t')
                qstart, qend, qseqid, sseqid = int(ls[2]), int(ls[3]), ls[0], ls[5]
                lineage = accession2lineage[sseqid.rsplit('_', 1)[0]] if 'plasmid@' not in sseqid else 'plasmid'
                if lineage == 'plasmid' and not plasmid:
                    continue

                ## filter out non-overlapping alignments
                if (AS := int(ls[14].split('AS:i:')[-1])) > scores[qseqid].get(lineage, -np.inf):
                    if any(compute_overlap((qstart, qend, *qcoord)) > 0 for qcoord in qcoords[qseqid]):
                        scores[qseqid][lineage] = AS
                        max_scores[qseqid] = max(max_scores.get(qseqid, -np.inf), AS)

                        alignments.append([qseqid, sseqid, AS, lineage])

        ## filter out low-score alignments
        duplicates = set()
        alignments.sort(key=lambda alignment: (alignment[0], alignment[2]))

        while alignments:
            alignment = alignments.pop()
            if max(alignment[2] / 0.995, alignment[2] + 50) > max_scores.get(alignment[0]):
                if (alignment[0], alignment[-1]) not in duplicates:
                    self.alignments.append(alignment)
                    duplicates.add((alignment[0], alignment[-1]))

    def run_mcl(self, max_iterations=1000, inflation=2, expansion=2):
        '''
        MCL for graph clustering.
        '''
        nodes = np.unique([hit[0] for hit in self.hits])
        node2index = {node: index for index, node in enumerate(nodes)}

        identities = {}
        while self.overlaps:
            overlap = self.overlaps.pop()
            if overlap[-1] <= self.DV:
                if (row := node2index.get(overlap[0])) and (col := node2index.get(overlap[1])):
                    if col < row:
                        row, col = col, row
                    identities[(row, col)] = max(1 - overlap[-1], identities.get((row, col), 0))

        rows, cols = zip(*list(identities.keys()))
        matrix = coo_matrix((list(identities.values()) * 2, (rows + cols, cols + rows)), shape=(len(nodes), len(nodes)))
        del self.overlaps, identities

        clusters = mcl(matrix, max_iterations=max_iterations, inflation=inflation, expansion=expansion)
        index2label = {index: label for label, indexes in enumerate(clusters) for index in indexes}
        labels = np.array([index2label.get(index) for index in range(len(nodes))])
        self.clusters = [nodes[labels == label] for label in np.unique(labels)[::-1]]

    def run_sc(self):
        '''
        Taxonomic assignments for each cluster with set cover.
        '''
        for cluster in self.clusters:
            elements = set(cluster)

            ## append scores of covered reads
            subsets = defaultdict(set)
            scores = defaultdict(lambda: defaultdict(dict))
            plasmid_qseqids, plasmid_lineages = set(), set()
            alignments = [alignment for alignment in self.alignments if alignment[0] in elements]
            for alignment in sorted(alignments, key=lambda alignment: (self.lineage2genome.get(alignment[-1], 0), alignment[-1]), reverse=True):
                if alignment[-1] == 'plasmid':
                    plasmid_qseqids.add(alignment[0])
                else:
                    subsets[alignment[-1]].add(alignment[0])
                    scores[alignment[-1]][alignment[0]] = alignment[2]

            ## get covers
            qseqid2lineage = {}
            if covers := set_cover(elements, subsets, scores):
                if len(covers) > 1:
                    score = {}
                    for qseqid in set.union(*[subsets[cover] for cover in covers]):
                        for cover in covers:
                            AS = scores[cover].get(qseqid, -np.inf)
                            if AS > score.get(qseqid, -np.inf):
                                qseqid2lineage[qseqid] = cover
                                score[qseqid] = AS
                else:
                    qseqid2lineage.update({qseqid: covers[0] for qseqid in subsets[covers[0]]})

            ## update for remaining elements
            qseqid2lineage.update({qseqid: 'unclassified' for qseqid in elements})

            ## check if over 50% reads of a species are labelled with plasmids
            for cover in set(qseqid2lineage.values()):
                qseqids = [qseqid in plasmid_qseqids for qseqid in {qseqid for qseqid, lineage in qseqid2lineage.items() if lineage == cover}]
                if sum(qseqids) / len(qseqids) > 0.5:
                    plasmid_lineages.add(cover)

            self.assignments.update({qseqid: (f'{lineage}@plasmid' if lineage in plasmid_lineages else lineage) for qseqid, lineage in qseqid2lineage.items()})

    def filter_hsp(self, subject_cover=90):
        '''
        Filter low-subject-cover and overlapping hits.
        '''
        scovs = defaultdict(lambda: defaultdict(set))
        bitscores, merged_bitscores = defaultdict(lambda: defaultdict(lambda: 0)), defaultdict(dict)
        for hit in self.hits:
            scovs[hit[0]][(hit[1], hit[5])].add((hit[6], hit[7]))
            bitscores[hit[0]][hit[1]] += hit[10]

        ## check whether reads collectively cover of a gene sequence
        discarded_hits = defaultdict(set)
        for elements in self.clusters:
            merged_scovs = defaultdict(set)
            for scov in [scovs.get(element) for element in elements]:
                for sseqid, coordinates in scov.items():
                    merged_scovs[sseqid].update(coordinates)

            discarded_sseqids = set()
            for sseqid, coordinates in merged_scovs.items():
                coordinates = list(sorted(coordinates))
                merged_coordinates = [list(coordinates[0])]

                for coordinate in coordinates[1:]:
                    if coordinate[0] <= merged_coordinates[-1][1]:
                        merged_coordinates[-1][1] = max(merged_coordinates[-1][1], coordinate[1])
                    else:
                        merged_coordinates.append(list(coordinate))

                ## discard hits they together < subject cover cutoff
                if sum(coordinate[1] - coordinate[0] for coordinate in merged_coordinates) / sseqid[1] < subject_cover / 100:
                    discarded_sseqids.add(sseqid[0])

            ## rank subtypes according to total bitscores
            sseqids = {sseqid[0] for sseqid in merged_scovs.keys()} - discarded_sseqids
            for sseqid in sseqids:
                subtype = sseqid.split('|')[2]
                bitscore = sum(bitscores.get(element).get(sseqid, 0) for element in elements)
                for element in elements:
                    merged_bitscores[element][subtype] = max(bitscore, merged_bitscores.get(element, {}).get(subtype, 0))

            for element in elements:
                discarded_hits[element].update(discarded_sseqids)

        nhits, self.hits = len(self.hits), [hit for hit in self.hits if hit[1] not in discarded_hits.get(hit[0], {})]

        ## filter out hits if they locate on the same position (>=25% overlap)
        hits = []
        qcoords = defaultdict(set)
        for hit in sorted([hit + [merged_bitscores.get(hit[0], {}).get(hit[1].split('|')[2], 0)] for hit in self.hits], key=lambda hit: (hit[-1], hit[-2]), reverse=True):
            qseqid, qstart, qend = hit[0], hit[3], hit[4]
            if (
                qseqid not in qcoords or
                all(compute_overlap((qstart, qend, *qcoord), max) < 0.25 for qcoord in qcoords[qseqid])
            ):
                qcoords[qseqid].add((qstart, qend))
                hits.append(hit)

        logger.info(
            f'... '
            f'read clusters: {len(self.clusters)} | '
            f'low-subject-cover HSPs: {nhits - len(self.hits)} | '
            f'overlapping HSPs: {len(self.hits) - len(hits)} | '
            f'remaining HSPs: {len(hits)}')

        self.hits = hits

    def run(self, debug=False, plasmid=False, skip_clean=False,
            max_target_seqs=25, evalue=1e-5, identity=0, subject_cover=90,
            secondary_num=2147483647, secondary_ratio=0.9,
            min_genome_copies=1, chunk_size=0,
            max_iterations=1000, inflation=2, expansion=2):
        '''
        Run the pipeline.
        '''
        self.debug = debug

        ## initial overlapping
        logger.info('Overlapping ...')
        self.run_overlap(mode='pre', identity=identity)

        ## annotating ARGs
        logger.info('Annotating ARGs ...')
        if not self.debug: self.run_diamond(max_target_seqs=max_target_seqs, evalue=evalue, identity=self.identity)
        self.parse_diamond()

        ## overlapping of ARG-containing reads
        logger.info('Overlapping of ARG-containing reads ...')
        self.run_overlap(mode='post', identity=identity, chunk_size=chunk_size)

        ## taxonomic assignment
        logger.info('Assigning taxonomy ...')
        if not self.debug: self.run_minimap(secondary_num=secondary_num, secondary_ratio=secondary_ratio)

        ## graph clustering
        logger.info('Graph clustering ...')
        self.run_mcl(max_iterations=max_iterations, inflation=inflation, expansion=expansion)
        self.filter_hsp(subject_cover=subject_cover)

        ## parse taxonomic assignments then reassign with set cover
        logger.info('Set covering ...')
        self.parse_minimap(plasmid=plasmid)
        self.run_sc()

        ## postprocessing
        reads = {hit[0]: {'remark': 'ARG-containing', 'hit': []} for hit in self.hits}
        lineage2copy = defaultdict(lambda: 0)
        for hit in self.hits:
            lineage = self.assignments.get(hit[0])
            carrier = 'plasmid' if '@plasmid' in lineage else 'chromosome'
            if self.lineage2genome.get(lineage.split('@')[0], -np.inf) < min_genome_copies:
                lineage = 'unclassified'
            else:
                lineage = lineage.split('@')[0]
            reads[hit[0]]['lineage'] = lineage
            reads[hit[0]]['plasmid'] = True if carrier == 'plasmid' else False

            reads[hit[0]]['hit'].append(hit[1].split('|', 1)[-1])
            lineage2copy[(lineage, *re.sub('@[A-Z-0-9]+', '', hit[1].replace('_', ' ')).split('|')[1:3], carrier)] += hit[9]

        self.profile = []
        for lineage, copy in lineage2copy.items():
            if (genome := self.lineage2genome.get(lineage[0], 0)) >= max(min_genome_copies, np.finfo(float).eps):
                self.profile.append([*lineage, copy, genome, copy / genome])
            else:
                self.profile.append([*lineage, copy, 0, 0])

        ## ARG profile
        with open(f'{self.outfile}.sarg.tsv', 'w') as w:
            if plasmid:
                w.write('\t'.join(['lineage', 'type', 'subtype', 'carrier', 'copy', 'genome', 'abundance']) + '\n')
                for row in sorted(self.profile):
                    w.write('\t'.join(row[:4] + [f'{row[4]:.3f}', f'{row[5]:.3f}', f'{row[6]:.3f}']) + '\n')
            else:
                w.write('\t'.join(['lineage', 'type', 'subtype', 'copy', 'genome', 'abundance']) + '\n')
                for row in sorted(self.profile):
                    w.write('\t'.join(row[:3] + [f'{row[4]:.3f}', f'{row[5]:.3f}', f'{row[6]:.3f}']) + '\n')

        ## read classification
        with open(f'{self.outfile}.sarg.json', 'w') as w:
            json.dump(dict(sorted(reads.items())), w, indent=4)

        ntypes = len({row[1] for row in self.profile if row[0] != 'unclassified'})
        nsubtypes = len({row[2] for row in self.profile if row[0] != 'unclassified'})
        nlineages = len({row[0] for row in self.profile if row[0] != 'unclassified'})
        abundance = sum(lineage2copy.values()) / sum(self.lineage2genome.values())

        logger.info(
            f'... ARG copies per genome: {abundance:.3f}'
            f' | ARG types: {ntypes}'
            f' | ARG subtypes: {nsubtypes}'
            f' | ARG-carrying species: {nlineages}.'
        )

        ## clean up
        if not skip_clean:
            for file in glob.glob(f'{self.outfile}.*.tmp'):
                os.remove(file)
