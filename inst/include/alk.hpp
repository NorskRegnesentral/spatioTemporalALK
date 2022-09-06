using namespace R_inla; //includes SPDE-spesific functions, e.g. Q_spde()
using namespace density; //use GMRF
using namespace Eigen; //Utilize sparse structures

template <class Type>
  Type nllALK(dataSet<Type> dat, paraSet<Type> par, spde_t<Type> spdeMatricesST_alk, LOSM_t<Type> A_alk_list){
    
    Type nll = 0;
    
    int nAges = dat.ageRange(1) - dat.ageRange(0) + 1;
    int minAge = dat.ageRange(0);
    int nObs = dat.length.size();
    
    Type sigma_beta0_alk = exp(par.log_sigma_beta0_alk(0));
    if(dat.rwBeta0_alk ==1){
      for(int y=1;y<dat.idx1.size();y++){
        for(int a=1; a<(nAges-1); ++a){
          nll -= dnorm(par.beta0_alk(y,a),par.beta0_alk(y-1,a-1),sigma_beta0_alk,true);
        }
      }
    }
    
    Type sigma =exp(par.logSigma_alk(0));
    Type kappa =exp(par.logKappa_alk(0));
    Type rho =2*invlogit(par.transRho_alk(0))-1;
    
    SparseMatrix<Type> Q = Q_spde(spdeMatricesST_alk,kappa);
    SparseMatrix<Type> Q_age(nAges,nAges);
    for(int a = 0; a< nAges; ++a){
      Q_age.coeffRef(a,a)=1;
    }
    
    nll += SEPARABLE(GMRF(Q_age),SEPARABLE(AR1(rho),GMRF(Q)))(par.xST_alk); //Opposite order than on R side
    
    Type scaleST = Type(1)/((4*3.14159265)*kappa*kappa); //No effect on results, but needed for interpreting the sigma^2 parameter as marginal variance. See section 2.1 in Lindgren (2011)
    
    matrix<Type> linPredMatrix(nObs, nAges-1);
    linPredMatrix.setZero(); 
    for(int a = 0; a<(nAges-1); ++a){
      for(int y =0; y<dat.idx1.size(); ++y){
        SparseMatrix<Type> A = A_alk_list(y);
        vector<Type> delta = (A*par.xST_alk.col(a).col(y).matrix())/sqrt(scaleST);
        
        linPredMatrix.col(a).segment(dat.idx1(y),dat.idx2(y)-dat.idx1(y)+1) = par.beta0_alk(y,a)+ 
            par.betaLength_alk(a)*dat.length.segment(dat.idx1(y),dat.idx2(y)-dat.idx1(y)+1) +
            delta*sigma;
        
//        linPredMatrix.col(a).segment(dat.idx1(y),dat.idx2(y)-dat.idx1(y)+1) = par.beta0_alk(y*(nAges-1) + a) + 
//          par.betaLength_alk(a)*dat.length.segment(dat.idx1(y),dat.idx2(y)-dat.idx1(y)+1) +
//          delta*sigma;
      } 
    }
    
    
    matrix<Type> ALK(nObs, nAges);
    ALK.setZero();
    
    vector<Type> tmp;
    vector<Type> probLess = ALK.col(0) ;
    for(int a = 0; a<nAges; ++a){
      probLess.setZero();
      for(int b = 0; b <a; ++b ){
        tmp = ALK.col(b);
        probLess = probLess + tmp;
      }
      if(a <(nAges-1)){
        tmp = linPredMatrix.col(a);
        ALK.col(a) = invlogit(tmp)*(1-probLess);
      }else{
        ALK.col(nAges-1) = (1-probLess);
      }
    }
    
    for(int s = 0; s<nObs; ++s){
      switch(dat.readability(s)){
      case 1:
        nll -= log(ALK(s,dat.age(s)- minAge));
        break;
      case 5:
        if(dat.age(s) <(nAges-1)){
          nll -= log(ALK(s,dat.age(s)- minAge) + ALK(s,dat.age(s)- minAge + 1) );
        }else{
          nll -= log(ALK(s,dat.age(s)- minAge));
        }
        break;
      case 6:
        if(dat.ageNotTruncated(s) < (nAges)){
          nll -= log(ALK(s,dat.age(s)- minAge) + ALK(s,dat.age(s)- minAge - 1) );
        }else{
          nll -= log(ALK(s,dat.age(s)- minAge));
        }
        break;
      }
    } 
    
    return(nll);
  }