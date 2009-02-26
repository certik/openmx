/* print used memory 

 2003.03.13 H. Kino
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define fp_bsize         1048576     /* buffer size for setvbuf */

static FILE *fp=NULL;
static char *filename=NULL;

typedef struct slist_{
  char *name;
  float size;
  struct slist_ *next;
} slist; 

static slist *top=NULL;



void PrintMemory(char *name, size_t size0, char *mode)
{
  slist *now, *last;
  float size;
  char buf[fp_bsize];          /* setvbuf */

  size=(float)size0;

/*  printf("PrintMemory called, %s %d %s\n",name,size,mode); */
  fflush(stdout);

  if (mode) {
    if (strcmp(mode,"init")==0) {

      /* try opening a file */
      if (filename) { free(filename); filename=NULL;}
      filename = strdup(name);
      fp=fopen(filename,"w");

#ifdef xt3
      setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

      if (fp==NULL) {
	printf("PrintMemory: fail to open a file: %s\n",filename);
	printf("PrintMemory: write to stdout\n"); 
	fp=stdout;
      }
      else {
	/* printf("PrintMemory: write to %s\n",filename); */
	fclose(fp);
      }


      return;
    }

    if (strcmp(mode,"sum")==0) {
      float total=0;
      for (now=top; now; now=now->next ) {
	total+= now->size;
      }

      if (fp==NULL) fp=stdout;
      if (fp!=stdout) {
	fp=fopen(filename,"a");

#ifdef xt3
        setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

	if (fp==NULL) {
	  fp=stdout;
	}
      }

      fprintf(fp,"Memory: %-35s %10.2f MBytes\n","total",total);

      if (fp!=stdout) {
	fclose(fp);
      }

      return;
    }

  }

  size /= 1024*1024;

  /* find iterator */
  for (last=now=top; now; now=now->next ) {
    if ( strcmp(now->name,name)==0 ) {
      /* found it*/
      break;
    }
    last=now;
  }
  if (now) { /* found */
    /* in this case, it is already written */
    return; 
  }
  else  { /* not found */
    if (last==NULL) { /* first time */
      last=top = (slist*) malloc(sizeof(slist));
    } 
    else { /* add this to the last */
      last=last->next=(slist*) malloc(sizeof(slist));
    }
    last->name=strdup(name);
    last->size=size;
    last->next=NULL;
  }

  /* write it to the file */

  if (fp==NULL) fp=stdout;
  if (fp!=stdout) {
    fp=fopen(filename,"a");

#ifdef xt3
    setvbuf(fp,buf,_IOFBF,fp_bsize);  /* setvbuf */
#endif

    if (fp==NULL) {
      fp=stdout;
    }
  }

  fprintf(fp,"Memory: %-35s %10.2f MBytes\n",name,size);

  if (fp!=stdout) {
    fclose(fp);
  }
 
}

#if 0

main()
{
/*   PrintMemory("test.memory",0,"init"); */
   PrintMemory("main: a",100000000,"add");
   PrintMemory("main: b",100000000,"add");
   PrintMemory("main: c",100000000,"add");
   PrintMemory("main: d",200000000,"add");
   PrintMemory("main: d",100000000,"add");
   PrintMemory("main: a",500000000,"add");
   PrintMemory("main: a",300000000,"add");
   PrintMemory("main: f",200000000,"add");

}
#endif







