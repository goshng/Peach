#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "sqlite3.h"
#include "kseq.h"

#define SIZE_SEQUENCE2 2
#define LOGLEVEL 0

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

  /* Create a list of strings.
   */
  int returnValue;
  sqlite3_stmt *stmt;
  sqlite3 *handle;

  /* Read each sequence name add add it to the table or the list.
   */
	secondFp = gzopen(argv[2], "r");
	secondSeq = kseq_init(secondFp);
  char sequenceName[1000];
  //int x1 = 0;
  int numberOfSequence2 = 0;
	int nextSequence2 = 0;
	while (nextSequence2 == 0) 
  {
    /* Create a new db at the start of a new section. */ 
    if (numberOfSequence2 == 0)
    {
      /* Read the second fastq file and add all of the names to the list 
       */
      returnValue = sqlite3_open(":memory:",&handle);
      char create_table[100] = "CREATE TABLE IF NOT EXISTS fastq "
                               "(name TEXT PRIMARY KEY,"
                               "comment TEXT)";
      returnValue = sqlite3_exec(handle,create_table,0,0,0);
      if(returnValue)
      {
        fprintf(stderr,"[%s] Info: Creating DB Failed\n", argv[0]);
        return 1;
      }
      else if (LOGLEVEL)
      {
        fprintf(stderr,"[%s] Info: Creating DB\n", argv[0]);
      }
    }

    /* Add second sequences. */
	  secondL = kseq_read(secondSeq);
	  if (secondL >= 0) {
      /* Insert sequences of FASTQ2. */
      numberOfSequence2++;
      sprintf(sequenceName, "INSERT INTO fastq VALUES('%s','%s')", 
              secondSeq->name.s, secondSeq->comment.s);
      returnValue = sqlite3_exec(handle,sequenceName,0,0,0);
      if(returnValue)
      {
        fprintf(stderr,"[%s] Error: Inserting data from DB Failed %s\n", 
                argv[0], secondSeq->name.s);
        return 1;
      }
      else if (LOGLEVEL)
      {
        fprintf(stderr,"[%s] Info: Inserting data from DB %s\n", 
                argv[0], secondSeq->name.s);
      }
    }
    else
    {
	    nextSequence2 = 1;
      if (LOGLEVEL)
      {
        fprintf(stderr,"[%s] Info: No more seq2\n", argv[0]);
      }
    }

    /* Change this number. */
    if (numberOfSequence2 == SIZE_SEQUENCE2 || nextSequence2 == 1)
    {
      /* Read the first fastq file and check if each name exists in the list.
       */
      firstFp = gzopen(argv[1], "r");
      firstSeq = kseq_init(firstFp);
      while ((firstL = kseq_read(firstSeq)) >= 0) {
        /* Check if this read exists in the list. */
        sprintf(sequenceName, "SELECT EXISTS(SELECT * FROM fastq WHERE "
                              "name='%s' LIMIT 1)",
                firstSeq->name.s);
        returnValue = sqlite3_prepare_v2(handle,sequenceName,-1,&stmt,0);
        if(returnValue)
        {
            fprintf(stderr, "Selecting data from DB Failed\n");
            return 1;
        }

        // Read the number of rows fetched
        int cols = sqlite3_column_count(stmt);
        //fprintf(stderr,"cols %d\n", cols); always 1
        returnValue = sqlite3_step(stmt);
            
        if(returnValue == SQLITE_ROW)
        {
          // Find the boolean of existence of the row.
          int col = 0;
          const char *val = (const char*)sqlite3_column_text(stmt,col);
          //printf("%s = %s\t",sqlite3_column_name(stmt,col),val);
          if (!strcmp(val,"1"))
          {
            printf("@%s %s\n", firstSeq->name.s, firstSeq->comment.s);
            printf("%s\n+\n%s\n", firstSeq->seq.s, firstSeq->qual.s);
            if (LOGLEVEL)
            {
              fprintf(stderr,"[%s] Info: selected sequence is %s\n",
                      argv[0], firstSeq->name.s);
            }
          }
//          else
//          {
//            fprintf(stderr,"[%s] Info: removed sequence is %s\n",
//                    argv[0], firstSeq->name.s);
//          }
        }
        else if(returnValue == SQLITE_DONE)
        {
          fprintf(stderr,"[%s] Error: impossible state.\n",
                  argv[0], firstSeq->name.s);
          return 1;
        }
        else
        {
          fprintf(stderr,"[%s] Error: unknown in sqlite.\n",argv[0]);
          return 1;
        }
      }
      kseq_destroy(firstSeq);
      gzclose(firstFp);

      /* Delete the list of strings.
       */
      sqlite3_close(handle);
      if (LOGLEVEL)
      {
        fprintf(stderr,"[%s] Info: deleting DB\n", argv[0]);
      }

      /* Reread the sequences in FASTQ2. */
      numberOfSequence2 = 0;
    }


	}
	kseq_destroy(secondSeq);
	gzclose(secondFp);


	return 0;
}
