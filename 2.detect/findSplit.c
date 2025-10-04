// findSplit.c -- find split sites in genome regions.
#include <stdio.h>
#include <stdlib.h>
#include "/home/songlizhi/software/htslib/htslib/sam.h"
int main(int argc, char * argv[])
{
	if (argc != 3) {
		printf("usage:    %s    bam    TEregion(chr:start-end)\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	samFile *fpin = sam_open(argv[1], "r");
	bam_hdr_t *hdr = sam_hdr_read(fpin);
	hts_idx_t *idx = sam_index_load(fpin, argv[1]);
	bam1_t *aln = bam_init1();
	hts_itr_t *itr = sam_itr_querys(idx, hdr, argv[2]);
	uint32_t *carray;
	int clen;

	if (fpin == NULL || hdr == NULL || idx == NULL || itr == NULL) {
		printf("Please check the path, header or index of your bamfile.");
		exit(EXIT_FAILURE);
	}

	while (sam_itr_next(fpin, itr, aln) >= 0) {
		carray = bam_get_cigar(aln);
		clen = aln->core.n_cigar;
		if (clen == 1) {
			continue;
		}
		else if (clen > 1) {
			//如果CIGAR第一位是softclip
			if (bam_cigar_opchr(carray[0]) == 'S') {
				printf("%s\t%ld\n", bam_get_qname(aln), aln->core.pos);
			}
			//如果CIGAR最后以为是softclip
			if (bam_cigar_opchr(carray[clen - 1]) == 'S') {
				printf("%s\t%ld\n", bam_get_qname(aln), aln->core.pos + bam_cigar2rlen(clen, carray));
			}
		}
	}

	bam_hdr_destroy(hdr);
	bam_itr_destroy(itr);
	hts_idx_destroy(idx);
	bam_destroy1(aln);
	sam_close(fpin);

	return 0;
}
