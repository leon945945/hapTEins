// getReadsTag.c -- get the alignment tag NM(Total number of mismatches and gaps) and AS(DP alignment score) of reads.
#include <stdio.h>
#include <stdlib.h>
#include "/home/songlizhi/software/htslib/include/htslib/sam.h"
#define bam_is_snd(b) (((b)->core.flag&BAM_FSECONDARY) != 0)
#define bam_is_sup(b) (((b)->core.flag&BAM_FSUPPLEMENTARY) != 0)
int main(int argc, char *argv[])
{
        if (argc != 3) {
                printf("usage:    %s    bamfile    region(chr:start-end)[.]\n", argv[0]);
                //exit(EXIT_FAILURE);
        }

        samFile *fp = sam_open(argv[1], "r");
        hts_idx_t *idx = sam_index_load(fp, argv[1]);
        bam_hdr_t *hdr = sam_hdr_read(fp);
        bam1_t *aln = bam_init1();
        hts_itr_t *itr = sam_itr_querys(idx, hdr, argv[2]);
        //uint32_t * cigar = 0;
        float ratio = 0.0;
        int alnLen = 0;
        int readLen = 0;
        char mark[4];
        //int i = 0;

        printf("readName\tchr\tbegin\tAS\tNM\talignLength\tcigarNum\treadLength\tratio\ttype\tMchr\tMbegin\n");
        while (sam_itr_next(fp, itr, aln) >= 0) {
                if (aln->core.tid < 0)
                        continue;
                //determine alignment type
                if (bam_is_sup(aln)) {
                        strcpy(mark, "sup");
                }
                else if (bam_is_snd(aln)) {
                        strcpy(mark, "snd");
                }
                else if (! bam_is_sup(aln) && ! bam_is_snd(aln)) {
                        strcpy(mark, "pri");
                }
                //get the alignment length and read length
                alnLen = bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln));
                readLen = aln->core.l_qseq;
                ratio = (float) alnLen/readLen;
                //output
                uint8_t *as_tag = bam_aux_get(aln, "AS");
                uint8_t *nm_tag = bam_aux_get(aln, "NM");
                printf("%s\t%s\t%ld\t%ld\t%ld\t%d\t%d\t%d\t%.2f\t%s\t%s\t%d",
                        bam_get_qname(aln),
                        hdr->target_name[aln->core.tid],
                        aln->core.pos,
                        as_tag ? bam_aux2i(as_tag) : -1,
                        nm_tag ? bam_aux2i(nm_tag) : -1,
                        alnLen,
                        aln->core.n_cigar,
                        readLen,
                        ratio,
                        mark,
                        hdr->target_name[aln->core.mtid],
                        aln->core.mpos);
                //cigar = bam_get_cigar(aln);
                //for (i = 0; i < aln->core.n_cigar; i++)
                        //printf("%c%d ", bam_cigar_opchr(cigar[i]), bam_cigar_oplen(cigar[i]));
                printf("\n");
        }

        bam_destroy1(aln);
        hts_itr_destroy(itr);
        bam_hdr_destroy(hdr);
        sam_close(fp);
}