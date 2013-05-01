#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	gzFile firstFp;
	kseq_t *firstSeq;
	gzFile secondFp;
	kseq_t *secondSeq;
	int firstL, secondL;
	if (argc == 2) {
		fprintf(stderr, "Usage: %s <r1.fastq.gz> <f2.fastq.gz>\n", argv[0]);
		return 1;
	}

  /* If the paired-end sequences are paired one by one,
     we could check just each sequence one by one and if the Illumina 
     filter is N for both pairs, then we can print them out. */ 
	firstFp = gzopen(argv[1], "r");
	firstSeq = kseq_init(firstFp);
	secondFp = gzopen(argv[2], "r");
	secondSeq = kseq_init(secondFp);
	while ((firstL = kseq_read(firstSeq)) >= 0) {
	  secondL = kseq_read(secondSeq);
    if (secondL < 0)
    {
      fprintf(stderr,"[%s] Error: No reads in the second file.\n",argv[0]);
      return 1;
    }

    /* Check if two reads share the same name. */
    if (strcmp(firstSeq->name.s,secondSeq->name.s))
    {
      fprintf(stderr,"[%s] Error: Two reads are not paired.\n",argv[0]);
      return 1;
    }

    if (firstSeq->comment.s[2] == 'N' && secondSeq->comment.s[2] == 'N') 
    {
		  printf("@%s %s\n", firstSeq->name.s, firstSeq->comment.s);
		  printf("%s\n+\n%s\n", firstSeq->seq.s, firstSeq->qual.s);
    }
	}
	secondL = kseq_read(secondSeq);
  if (secondL >= 0)
  {
    fprintf(stderr,"[%s] Error: More reads in the second file.\n",argv[0]);
    return 1;
  }

	kseq_destroy(firstSeq);
	gzclose(firstFp);
	kseq_destroy(secondSeq);
	gzclose(secondFp);

	return 0;
}
