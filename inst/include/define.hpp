
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



#define ADREPORT_F(name,F) F->reportvector.push(name,#name);

#define SIMULATE_F(F)				\
if(isDouble<Type>::value && F->do_simulate)


#define REPORT_F(name,F)					                               \
if(isDouble<Type>::value && F->current_parallel_region<0) {	\
  Rf_defineVar(Rf_install(#name),					                      \
               PROTECT(asSEXP(name)),F->report);			         \
  UNPROTECT(1);						                                       \
}




template <class Type>
  struct dataSet;

template <class Type>
  struct paraSet;


template <class Type>
struct dataSet{
  vector<Type> age; 
  vector<int> ageNotTruncated; 
  vector<Type> length; 
  vector<int> readability; 
  vector<int> ageRange; 
  vector<int> idx1; 
  vector<int> idx2; 
  int rwBeta0_alk;
  int maxAge;
  int minAge;
  int usePCpriorsALK;
  int spatioTemporalALK;
  int spatialALK;
  int betaLength;
  vector<Type> pcPriorsALKRange; 
  vector<Type> pcPriorsALKSD; 
//  LOSM_t<Type> A_alk_list;

};


template <class Type>
struct paraSet{
  matrix<Type> beta0_alk; 
  vector<Type> log_sigma_beta_alk;
  vector<Type> betaLength_alk; 
  vector<Type> logSigma_alk; 
  vector<Type> logKappa_alk; 
  vector<Type> transRho_alk; 
  array<Type> xST_alk; 
  array<Type> xS_alk; 
  
};
