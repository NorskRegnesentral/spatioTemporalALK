
/* List of sparse matrices */
template<class Type>
struct LOSM_t : vector<SparseMatrix<Type> > {
  LOSM_t(SEXP x){  /* x = List passed from R */
(*this).resize(LENGTH(x));
    for(int i=0; i<LENGTH(x); i++){
      SEXP sm = VECTOR_ELT(x, i);
      if(!isValidSparseMatrix(sm))
        error("Not a sparse matrix");
      (*this)(i) = asSparseMatrix<Type>(sm);
    }
  }
};


template <class Type>
  struct dataSet;

template <class Type>
  struct paraSet;


template <class Type>
struct dataSet{
  vector<int> age; 
  vector<int> ageNotTruncated; 
  vector<Type> length; 
  vector<int> readability; 
  vector<int> ageRange; 
  vector<int> idx1; 
  vector<int> idx2; 
  int rwBeta0_alk;
//  LOSM_t<Type> A_alk_list;

};


template <class Type>
struct paraSet{
  matrix<Type> beta0_alk; 
  vector<Type> log_sigma_beta0_alk;
  vector<Type> betaLength_alk; 
  vector<Type> logSigma_alk; 
  vector<Type> logKappa_alk; 
  vector<Type> transRho_alk; 
  array<Type> xST_alk; 
  
};
