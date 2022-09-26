#define TMB_LIB_INIT R_init_spatioTemporalALK

#include <TMB.hpp>
using namespace tmbutils;
#include "../inst/include/define.hpp"
#include "../inst/include/alk.hpp"
 


template<class Type>  
  Type objective_function<Type>::operator() ()
{
  using namespace density; //use GMRF
  using namespace Eigen; //Utilize sparse structures
  using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde() 
  
  dataSet<Type> dataset;  
  DATA_IVECTOR(age); dataset.age = age;//The response
  DATA_IVECTOR(ageNotTruncated); dataset.ageNotTruncated = ageNotTruncated;
  DATA_VECTOR(length); dataset.length = length;//Covariate
  DATA_IVECTOR(readability); dataset.readability = readability;
  DATA_IVECTOR(ageRange);dataset.ageRange = ageRange;  
  DATA_IVECTOR(idx1); dataset.idx1 = idx1;//Bookkeeping
  DATA_IVECTOR(idx2); dataset.idx2 = idx2;
  DATA_STRUCT(spdeMatricesST_alk,spde_t); //TODO: Include in dataset
  DATA_STRUCT(A_alk_list, LOSM_t); 
  DATA_STRUCT(A_alk_obs, LOSM_t); 
  DATA_INTEGER(rwBeta0_alk); dataset.rwBeta0_alk = rwBeta0_alk;
  DATA_INTEGER(maxAge); dataset.maxAge = maxAge;
  DATA_INTEGER(minAge); dataset.minAge = minAge;
  DATA_VECTOR(pcPriorsALKRange); dataset.pcPriorsALKRange = pcPriorsALKRange;
  DATA_VECTOR(pcPriorsALKSD); dataset.pcPriorsALKSD = pcPriorsALKSD;
  DATA_INTEGER(usePCpriorsALK); dataset.usePCpriorsALK = usePCpriorsALK;
  DATA_INTEGER(spatioTemporalALK); dataset.spatioTemporalALK = spatioTemporalALK;
  DATA_INTEGER(spatialALK); dataset.spatialALK = spatialALK;
  
  paraSet<Type> paraset; 
  PARAMETER_MATRIX(beta0_alk); paraset.beta0_alk = beta0_alk;//Intercepts
  PARAMETER_VECTOR(log_sigma_beta0_alk);paraset.log_sigma_beta0_alk = log_sigma_beta0_alk;      
  PARAMETER_VECTOR(betaLength_alk); paraset.betaLength_alk = betaLength_alk;//Regression parameters
  PARAMETER_VECTOR(logSigma_alk);paraset.logSigma_alk = logSigma_alk;
  PARAMETER_VECTOR(logKappa_alk);paraset.logKappa_alk = logKappa_alk;
  PARAMETER_VECTOR(transRho_alk);paraset.transRho_alk = transRho_alk;
  PARAMETER_ARRAY(xS_alk); paraset.xS_alk = xS_alk;// Ordering: xST.col(age).col(year).col(spatial)
  PARAMETER_ARRAY(xST_alk); paraset.xST_alk = xST_alk;// Ordering: xST.col(age).col(year).col(spatial)
  
  
  Type nll = 0;
    
  nll += nllALK(dataset,paraset,spdeMatricesST_alk,A_alk_list); 
  
  
  return nll;
}
