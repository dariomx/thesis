#include <cs.h>


void sccblocks_(int *row, int *col, double *val, 
		int *nPtr, int *nelemsPtr, 
		int *nbPtr, int *p, int *r, 
		int* brows, int* bcols, double* bvals,
		int* brangelimits, int* bnnzlimits) {
  /*
    row
    INPUT Array of integer row indices i for matrix A
    Size: nelems.
    
    col
    INPUT Array of integer column indices j, for matrix A. 
    Size: nelems.
    
    val
    INPUT Array of corresponding double precision values A(i,j). 
    Size: nelems.
	
    nPtr
    INPUT Pointer to the integer size n of each dimension in (square) matrix A.
    An integer in "Fortran side".
    
    nelemsPtr
    INPUT Pointer to the integer number of elements nelems in matrix A.
    An integer in "Fortran side".
    
    nbPtr
    OUTPUT Pointer to the integer number nb of strongly connected components (number of blocks).
    An integer in "Fortran side".
    
    p 
    OUTPUT Array of integers; this is the permutation vector for revealing the blocks;
    the permuted matrix Ap reads as Ap = A(p, p) in Matlab notation. 
    Size: n.
    
    r 
    OUTPUT Array of integers; pairs of successive elements in r  denote the [start, end) limits 
    of the blocks in the permuted matrix Ap. 
    Size: n + 1; however only nb + 1 entries will be valid on returning (nb <= n).
    
    brows
    OUTPUT Array of integers; row indices i for all generated blocks, 
    for the consecutive blocks in Ap, collected in a single array. 
    Size: nelems; however only "total number of nonzeros in generated blocks" entries will be valid on returning.
     
    bcols
    OUTPUT Array of integers; column indices j for all generated blocks, 
    for the consecutive blocks in Ap, collected in a single array. 
    Size: nelems; however only "total number of nonzeros in generated blocks" entries will be valid on returning.
    
    bvals
    OUTPUT Array of  double precision values at corresponding (i, j) locations for all generated blocks,
    collected in a single array. 
    Size:  nelems; however only "total number of nonzeros in generated blocks" entries will be valid on returning. 
    
    brangelimits
    OUTPUT Array of integer values; the range limits of each block in the permuted matrix Ap.
    The ith block will be located in the [ brangelimits[i], brangelimits[i+1] ) range for
    both row and column indices in Ap. 
    It follows that the size of the ith block can be computed by subtracting:
    brangelimits[i+1] - brangelimits[i].
    Size: n + 1; however only nb + 1 entries will be valid on returning (nb <= n).
     
    bnnzlimits
    OUTPUT Array of integer values; consecutive entries accumulate the nnz's for the successive blocks in Ap.
    This is useful for finding the index ranges for traversals in brows[] 
    (and similarly for bcols[] and bvals[]) for each block.
    So for the ith block we should use indices k in range [ bnnzlimits[i], bnnzlimits[i+1] )
    for accessing brows[k], bcols[k], bvals[k] CO format triplets.
    Note that indices in brows[], bcols[] for each block matrix start at 1;
    so each block matrix is NOT assumed as part of Ap, it can be readily extracted and used
    independently.
    Size: n + 1; however only nb + 1 entries will be valid on returning (nb <= n);
  */


  
  /* Matrix related variables */
  cs *triplets, *A, *Ap;
  int n = *nPtr;
  int nelems = *nelemsPtr;
  int nnz;

  /* Strongly connected components (scc) related variables  */
  csd *sccinfo; 
  int nb;
  csi csinb;
  csi *csip, *csir, *csipinv;
    
  /* Block-generation related variables*/
  cs **tblocks, **cblocks;
  int *bsizes, *bnnzs;
  
  int i;
  int ib;
  csi j, k;
  csi bstart, bend, bsize;
  int bnnz, bcol;
  csi loc;
    
    
  /* Build the graph, taking care of the the fact that row[], col[] 
     come from Fortran so they assume 1-based numbering (input arguments from the "Fortran side") 
  */
  /* XXX casting */
  triplets = cs_spalloc(n, n, nelems, (csi)1, (csi)1);
  for(i = 0; i < nelems; i++) {
    cs_entry(triplets, row[i]-1, col[i]-1, val[i]);
  }
  
  
  /* Convert the matrix for COordinate format (CO, i.e. triplets) to Compressed Sparse Column (CSC) format.
     This is mandatory in CSparse in order to have access to its full set of
     matrix operations.
  */
  A = cs_compress(triplets);
  
  
  /* Remove zero elements; remove duplicates.
     Note that Matrix Market format read facilities provided for C and Fortran do not
     automatically drop zero values from the triplets of entries in a "coordinate format" file.
     However Matlab's mmread() facility automatically drops these zero values because it
     uses sparse() to generate the final matrix; sparse() actually drops these zeros.
     So nnz is the number of actual non-zero elements in A.
  */
  cs_dropzeros(A);
  cs_dupl(A);
  nnz = 0;
  for(j=0; j < n; j++) {
    loc = A->p[j];
    for( ; loc < A->p[j+1]; loc++) {
      nnz++;
    }
  }
  
  /* Compute strongly connected components 
   */
  sccinfo = cs_scc(A); 
  
  
  /* Extract the relevant info from the returned structure:
     Also taking care of the 0- to 1- based numbering for p[] , r[] since these
     will be used in the Fortran side (output arguments to the "Fortran side").
  */
  csinb = sccinfo->nb;
  /* XXX casting */
  nb = (int)csinb;
  *nbPtr = nb;
    
  csip = (csi*) malloc(n * sizeof(csi));
  csir = (csi*) malloc((nb + 1) * sizeof(csi));
  for(i = 0; i < n; i++) {
    /* XXX casting */
    p[i] = (int)sccinfo->p[i] + 1;
    csip[i] = sccinfo->p[i];
  }
    
  for(i = 0; i < nb + 1; i++) {
    /* XXX casting */
    r[i] = (int)sccinfo->r[i] + 1;
    csir[i] = sccinfo->r[i];
  }
  
  
  /* In order to ask CSparse to apply the permutation and get Ap = A(p, p)
     we also need the inverse permutation pinv. 
     We use the name csipinv for this.
     We generally prepend the names of other variables as well with csi.
     csi in CSparse is used as an alias for the ptrdiff_t data type.
     ptrdiff_t is a signed integer type which holds the difference between
     two pointers. 
     XXX: Currently we cast ptrdiff_t to int in various places but we should look more carefully into this!
     sizeof(ptrdiff_t)= 8 bytes while sizeof(int)= 4 bytes in a couple of platforms I checked;
     this means we must be careful when casting csi <-> int in the following cases:
     i)  if outside the range [-2^(31), 2^(31)-1] i.e. a couple of billions.
     ii) in operations involving sizeof, e.g. malloc()
     Look for the lines prepended with the "XXX casting" comment.
  */
  csipinv = (csi*) malloc(n * sizeof(csi));
  for(i=0; i < n; i++) {
    /* XXX casting */
    csipinv[csip[i]] = (csi)i;
  }

  
  /* Compute Ap = A(p, p) */
  /* XXX casting */
  Ap = cs_permute(A, csipinv, csip, (csi)1);
  
  
  /* Allocate space for the list of matrix block representations: 
     tblocks[][] : a list of CO format representations (in CSparse)
     cblocks[][] : a list of CSC representations (in CSparse)
  */   
  tblocks = (cs**) malloc(nb * sizeof(cs*));
  cblocks = (cs**) malloc(nb * sizeof(cs*));

  
  /* Allocate space for:
     bsizes[] : the array of matrix block sizes
     bnnzs[]   : the array of nnz for the matrix blocks
     One array entry per consecutive block in Ap
  */
  bsizes = (int*) malloc(nb * sizeof(int));
  bnnzs = (int*) malloc(nb * sizeof(int));
  
  
  /* Loop for generating the blocks */
  for(ib=0; ib < csinb; ib++) {
    bstart = csir[ib];
    bend = csir[ib+1];
    bsize = bend - bstart;
    bsizes[ib] = bsize;
    bnnz = 0;
    /* First pass: loop for counting non-zeros in the matrix block;
       store them in  bnnzs[] in the process. */
    for(j=bstart; j < bend; j++) {
      loc = Ap->p[j];
      for( ; loc < Ap->p[j+1]; loc++) {
	if((Ap->i[loc] >= bstart) && (Ap->i[loc] < bend)) {
	  bnnz++;
	} 
      }
    }
    bnnzs[ib] = bnnz;
    
    /* Allocate the block */
    tblocks[ib] = cs_spalloc(bsize, bsize, bnnz, (csi)1, (csi)1);
    
    /* Second pass: loop for filling the block. */
    for(j=bstart, bcol=0; j < bend; j++) {
      loc = Ap->p[j];
      for( ; loc < Ap->p[j+1]; loc++) {
	if((Ap->i[loc] >= bstart) && (Ap->i[loc] < bend)) {
	  cs_entry(tblocks[ib], Ap->i[loc] - bstart, bcol, Ap->x[loc]);
	} 
      }
      bcol++;
    }
    cblocks[ib] = cs_compress(tblocks[ib]);
  }
  
  
  /* Fill in remaining block-related output arguments
     Also take care of possible 0- to 1-based issues in the process 
  */
  i = 0;
  for(ib=0; ib < csinb; ib++) {
    for(k = 0; k < tblocks[ib]->nz; k++) {
      brows[i] = (tblocks[ib]->i[k]) + 1;
      bcols[i] = (tblocks[ib]->p[k]) + 1;
      bvals[i] = tblocks[ib]->x[k];
      i++;
    }
  }
  
  brangelimits[0] = 1;
  for(ib=1; ib <= nb; ib++) {
    brangelimits[ib] = brangelimits[ib-1] + bsizes[ib-1];
  }  
  
  bnnzlimits[0] = 1;
  for(ib=1; ib <= nb; ib++) {
    bnnzlimits[ib] = bnnzlimits[ib-1] + bnnzs[ib-1];
  }
  

  /* free dynamically allocated space */
  cs_spfree(triplets);
  for(ib=0; i< nb; i++) {
    cs_spfree(tblocks[ib]);
    cs_spfree(cblocks[ib]);
  }
  cs_spfree(A);
  cs_spfree(Ap);
  free(csip);
  free(csir);
  free(csipinv);
  free(tblocks);
  free(cblocks);
  free(bsizes);
  free(bnnzs);
    
}
