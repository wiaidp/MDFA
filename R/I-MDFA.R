#' Simplified call for multivariate MSE estimation: no constraints, no regularization, no customization
#'
#' @param L Filter-length
#' @param weight_func Spectrum: DFT-matrix or alternative (for example model-based) estimate: first column is the target variable, additional columns are explanatory variables
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param Gamma Generic target specification: typically symmetric Lowpass (trend) or Bandpass (cycle) filters. Highpass and anticipative allpass (forecast) can be specified too
#' @return mdfa_obj MDFA object
#' @export
#'


MDFA_mse<-function(L,weight_func,Lag,Gamma)
{

  cutoff<-pi
  lin_eta<-F
  lambda<-0
  eta<-0
  weight_constraint<-rep(1/(ncol(weight_func)-1),ncol(weight_func)-1)
  lambda_cross<-lambda_smooth<-0
  lambda_decay<-c(0,0)
  lin_expweight<-F
  shift_constraint<-rep(0,ncol(weight_func)-1)
  grand_mean<-F
  b0_H0<-NULL
  c_eta<-F
  weights_only<-F
  weight_structure<-c(0,0)
  white_noise<-F
  synchronicity<-F
  lag_mat<-matrix(rep(0:(L-1),ncol(weight_func)),nrow=L)
  troikaner<-F
  i1<-i2<-F

  mdfa_obj<-mdfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                          lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                          b0_H0,c_eta,weight_structure,white_noise,
                          synchronicity,lag_mat,troikaner)

  return(list(mdfa_obj=mdfa_obj))
}





#' Simplified call for multivariate MSE estimation with filter constraints: no regularization, no customization
#'
#' @param L Filter-length
#' @param weight_func Spectrum: DFT-matrix or alternative (for example model-based) estimate: first column is the target variable, additional columns are explanatory variables
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param Gamma Generic target specification: typically symmetric Lowpass (trend) or Bandpass (cycle) filters. Highpass and anticipative allpass (forecast) can be specified too
#' @param i1 Boolean. If T a first-order filter constraint in frequency zero is obtained: amplitude of real-time filter must match weight_constraint (handles integration order one)
#' @param i2 Boolean. If T a second-order filter constraint in frequency zero is obtained: time-shift of real-time filter must match target (together with i1 handles integration order two)
#' @param weight_constraint Vector of amplitude constraints in frequency zero (typically 1 for lowpass and zero for bandpass or highpass)
#' @param shift_constraint Vector of time-shift constraints in frequency zero (typically zero for lowpass)
#' @return mdfa_obj MDFA object
#' @export
#'

MDFA_mse_constraint<-function(L,weight_func,Lag,Gamma,i1,i2,weight_constraint,shift_constraint)
{

  cutoff<-pi
  lin_eta<-F
  lambda<-0
  eta<-0
  lambda_cross<-lambda_smooth<-0
  lambda_decay<-c(0,0)
  lin_expweight<-F
  grand_mean<-F
  b0_H0<-NULL
  c_eta<-F
  weights_only<-F
  weight_structure<-c(0,0)
  white_noise<-F
  synchronicity<-F
  lag_mat<-matrix(rep(0:(L-1),ncol(weight_func)),nrow=L)
  troikaner<-F

  mdfa_obj<-mdfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                          lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                          b0_H0,c_eta,weight_structure,white_noise,
                          synchronicity,lag_mat,troikaner)

  return(list(mdfa_obj=mdfa_obj))
}





#' Simplified call for multivariate customization: no regularization, no constraints
#'
#' @param L Filter-length
#' @param weight_func DFT-matrix or alternative (for example model-based) estimate: first column is the target variable, additional columns are explanatory variables
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param Gamma Generic target specification: typically symmetric Lowpass (trend) or Bandpass (cycle) filters. Highpass and anticipative allpass (forecast) can be specified too
#' @param cutoff Specifies start-frequency in stopband from which Smoothness is emphasized (corresponds typically to the cutoff of the lowpass target). Is not used if eta=0.
#' @param lambda Customization parameter: Timeliness is emphasized in the ATS-trilemma if lambda>0
#' @param eta Customization parameter: Smoothness is emphasized in the ATS-trilemma if eta>0
#' @return mdfa_obj MDFA object
#' @export
#'

MDFA_cust<-function(L,weight_func,Lag,Gamma,cutoff,lambda,eta)
{

  lin_eta<-F
  weight_constraint<-rep(1/(ncol(weight_func)-1),ncol(weight_func)-1)
  lambda_cross<-lambda_smooth<-0
  lambda_decay<-c(0,0)
  lin_expweight<-F
  shift_constraint<-rep(0,ncol(weight_func)-1)
  grand_mean<-F
  b0_H0<-NULL
  c_eta<-F
  weights_only<-F
  weight_structure<-c(0,0)
  white_noise<-F
  synchronicity<-F
  lag_mat<-matrix(rep(0:(L-1),ncol(weight_func)),nrow=L)
  troikaner<-F
  i1<-i2<-F

  mdfa_obj<-mdfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                          lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                          b0_H0,c_eta,weight_structure,white_noise,
                          synchronicity,lag_mat,troikaner)

  return(list(mdfa_obj=mdfa_obj))
}







#' Simplified call for multivariate customization with filter constraints: no regularization
#'
#' @param L Filter-length
#' @param weight_func DFT-matrix or alternative (for example model-based) estimate: first column is the target variable, additional columns are explanatory variables
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param Gamma Generic target specification: typically symmetric Lowpass (trend) or Bandpass (cycle) filters. Highpass and anticipative allpass (forecast) can be specified too
#' @param cutoff Specifies start-frequency in stopband from which Smoothness is emphasized (corresponds typically to the cutoff of the lowpass target). Is not used if eta=0.
#' @param i1 Boolean. If T a first-order filter constraint in frequency zero is obtained: amplitude of real-time filter must match weight_constraint (handles integration order one)
#' @param i2 Boolean. If T a second-order filter constraint in frequency zero is obtained: time-shift of real-time filter must match target (together with i1 handles integration order two)
#' @param weight_constraint Vector of amplitude constraints in frequency zero (typically 1 for lowpass and zero for bandpass or highpass)
#' @param shift_constraint Vector of time-shift constraints in frequency zero (typically zero for lowpass)
#' @param lambda Customization parameter: Timeliness is emphasized in the ATS-trilemma if lambda>0
#' @param eta Customization parameter: Smoothness is emphasized in the ATS-trilemma if eta>0
#' @return mdfa_obj MDFA object
#' @export
#'

MDFA_cust_constraint<-function(L,weight_func,Lag,Gamma,cutoff,lambda,eta,i1,i2,weight_constraint,shift_constraint)
{

  lin_eta<-F
  lambda_cross<-lambda_smooth<-0
  lambda_decay<-c(0,0)
  lin_expweight<-F
  grand_mean<-F
  b0_H0<-NULL
  c_eta<-F
  weights_only<-F
  weight_structure<-c(0,0)
  white_noise<-F
  synchronicity<-F
  lag_mat<-matrix(rep(0:(L-1),ncol(weight_func)),nrow=L)
  troikaner<-F

  mdfa_obj<-mdfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                          lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                          b0_H0,c_eta,weight_structure,white_noise,
                          synchronicity,lag_mat,troikaner)

  return(list(mdfa_obj=mdfa_obj))
}









#' Simplified call for multivariate MSE and customization with regularization but no constraints
#'
#' @param L Filter-length
#' @param weight_func DFT-matrix or alternative (for example model-based) estimate: first column is the target variable, additional columns are explanatory variables
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param Gamma Generic target specification: typically symmetric Lowpass (trend) or Bandpass (cycle) filters. Highpass and anticipative allpass (forecast) can be specified too
#' @param cutoff Specifies start-frequency in stopband from which Smoothness is emphasized (corresponds typically to the cutoff of the lowpass target). Is not used if eta=0.
#' @param lambda Customization parameter: Timeliness is emphasized in the ATS-trilemma if lambda>0
#' @param eta Customization parameter: Smoothness is emphasized in the ATS-trilemma if eta>0
#' @param lambda_cross Regularization: cross-sectional term
#' @param lambda_decay Regularization: decay term
#' @param lambda_smooth Regularization: smoothness term
#' @param troikaner Boolean: if T then additional statistics such as ATS and degrees of freedom are computed (numerically involving)
#' @return mdfa_obj MDFA object
#' @export
#'

MDFA_reg<-function(L,weight_func,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth,troikaner=F)
{


  lin_eta<-F
  weight_constraint<-rep(1/(ncol(weight_func)-1),ncol(weight_func)-1)
  lin_expweight<-F
  shift_constraint<-rep(0,ncol(weight_func)-1)
  grand_mean<-F
  b0_H0<-NULL
  c_eta<-F
  weights_only<-F
  weight_structure<-c(0,0)
  white_noise<-F
  synchronicity<-F
  lag_mat<-matrix(rep(0:(L-1),ncol(weight_func)),nrow=L)
  i1<-i2<-F



  mdfa_obj<-mdfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                          lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                          b0_H0,c_eta,weight_structure,white_noise,
                          synchronicity,lag_mat,troikaner)

  return(list(mdfa_obj=mdfa_obj))
}









#' Simplified call for multivariate MSE/customization with regularization and constraints
#'
#' @param L Filter-length
#' @param weight_func DFT-matrix or alternative (for example model-based) estimate: first column is the target variable, additional columns are explanatory variables
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param Gamma Generic target specification: typically symmetric Lowpass (trend) or Bandpass (cycle) filters. Highpass and anticipative allpass (forecast) can be specified too
#' @param cutoff Specifies start-frequency in stopband from which Smoothness is emphasized (corresponds typically to the cutoff of the lowpass target). Is not used if eta=0.
#' @param i1 Boolean. If T a first-order filter constraint in frequency zero is obtained: amplitude of real-time filter must match weight_constraint (handles integration order one)
#' @param i2 Boolean. If T a second-order filter constraint in frequency zero is obtained: time-shift of real-time filter must match target (together with i1 handles integration order two)
#' @param lambda Customization parameter: Timeliness is emphasized in the ATS-trilemma if lambda>0
#' @param eta Customization parameter: Smoothness is emphasized in the ATS-trilemma if eta>0
#' @param weight_constraint Vector of amplitude constraints in frequency zero (typically 1 for lowpass and zero for bandpass or highpass)
#' @param shift_constraint Vector of time-shift constraints in frequency zero (typically zero for lowpass)
#' @param lambda_cross Regularization: cross-sectional term
#' @param lambda_decay Regularization: decay term
#' @param lambda_smooth Regularization: smoothness term
#' @param troikaner Boolean: if T then additional statistics such as ATS and degrees of freedom are computed (numerically involving)
#' @return mdfa_obj MDFA object
#' @export
#'

MDFA_reg_constraint<-function(L,weight_func,Lag,Gamma,cutoff,lambda,eta,lambda_cross,lambda_decay,lambda_smooth,i1,i2,weight_constraint,shift_constraint,troikaner=F)
{

  lin_eta<-F
  lin_expweight<-F
  grand_mean<-F
  b0_H0<-NULL
  c_eta<-F
  weights_only<-F
  weight_structure<-c(0,0)
  white_noise<-F
  synchronicity<-F
  lag_mat<-matrix(rep(0:(L-1),ncol(weight_func)),nrow=L)


  mdfa_obj<-mdfa_analytic(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                          lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                          b0_H0,c_eta,weight_structure,white_noise,
                          synchronicity,lag_mat,troikaner)

  return(list(mdfa_obj=mdfa_obj))
}






#' Main estimation routine: Sets-up the generic optimization criteria proposed in MDFA-Legacy project (book)
#'
#' @param L Filter-length
#' @param lambda Customization parameter: Timeliness is emphasized in the ATS-trilemma if lambda>0
#' @param weight_func DFT-matrix or alternative (for example model-based) estimate: first column is the target variable, additional columns are explanatory variables
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param Gamma Generic target specification: typically symmetric Lowpass (trend) or Bandpass (cycle) filters. Highpass and anticipative allpass (forecast) can be specified too
#' @param eta Customization parameter: Smoothness is emphasized in the ATS-trilemma if eta>0
#' @param cutoff Specifies start-frequency in stopband from which Smoothness is emphasized (corresponds typically to the cutoff of the lowpass target). Is not used if eta=0.
#' @param i1 Boolean. If T a first-order filter constraint in frequency zero is obtained: amplitude of real-time filter must match weight_constraint (handles integration order one)
#' @param i2 Boolean. If T a second-order filter constraint in frequency zero is obtained: time-shift of real-time filter must match target (together with i1 handles integration order two)
#' @param weight_constraint Constraint vector in the case i1==T
#' @param lambda_cross Regularization: cross-sectional term
#' @param lambda_decay Regularization: decay term
#' @param lambda_smooth Regularization: smoothness term
#' @param lin_eta Boolean: impose continuous or discontinuous Smoothness customization
#' @param shift_constraint Constraint vector in the case i2==T
#' @param grand_mean Boolean: if T then a grand-mean parametrization is imposed (default is F)
#' @param b0_H0 Regularization: shrinkage target (arbitrary designs can be replicated by imposing strong regularization)
#' @param c_eta Boolean: impose mild/strong smoothness customization (default is F)
#' @param weight_structure Add structure to the optimization criterion (default value is weight_structure<-c(0,0): no structure imposed)
#' @param white_noise Impose a flat DFT (phase information is maintained)
#' @param synchronicity Impose a zero shift across series (amplitude information is maintained)
#' @param lag_mat Matrix for implementing effective lags in a mixed-frequency setting
#' @param troikaner Boolean: if T then degrees oif freedom will be computed (time-consuming computations if K is large)
#' @return b Matrix of optimal filter coefficients
#' @return trffkt Complex transfer function of optimal multivariate filter
#' @return rever Criterion value (corresponds to a sample estimate of the MSE if lambda=eta=0)
#' @return degrees_freedom Degrees of freedom (when imposing regularization the degrres of freedom are smaller than L times the number of explanatory series)
#' @return Accuracy Accuracy term in decomposition of MSE
#' @return Smoothness Smoothness term in decomposition of MSE
#' @return Timeliness Timeliness term in decomposition of MSE
#' @return MS_error Sample estimate of MSE: consistent estimate in the cases lambda>0 and/or eta>0
#' @return freezed_degrees_new The complementary of degrees_freedom
#' @return tic Troikaner information criterion: based on sub-dimension of shrinkage space
#' @export
#'

mdfa_analytic<-function(L,lambda,weight_func,Lag,Gamma,eta,cutoff,i1,i2,weight_constraint,
                        lambda_cross,lambda_decay,lambda_smooth,lin_eta,shift_constraint,grand_mean,
                        b0_H0,c_eta,weight_structure,white_noise,
                        synchronicity,lag_mat,troikaner)
{

# Enforce meaningful parameter values
  lambda<-abs(lambda)
  eta<-abs(eta)
  cutoff<-min(abs(cutoff),pi)

# Resolution of discrete frequency-grid
  K<-nrow(weight_func)-1
  weight_target<-weight_func[,1]
  weight_func<-weight_func*exp(-1.i*Arg(weight_target))

  white_noise_synchronicity_obj<-white_noise_synchronicity(weight_func,white_noise,synchronicity)

  weight_func<-white_noise_synchronicity_obj$weight_func

  if (i2&i1&any(weight_constraint==0))
  {
    print(rep("!",100))
    print("Warning: i2<-T is not meaningful when i1<-T and weight_constraint=0 (bandpass)")
    print(rep("!",100))
  }
  if (!(length(b0_H0)==L*(dim(weight_func)[2]-1))&length(b0_H0)>0)
    print(paste("length of b0_H0 vector is ",length(b0_H0),": it should be ",L*(dim(weight_func)[2]-1)," instead",sep=""))
# The function spect_mat_comp rotates all DFTs such that first column is real!#Lag<-0

  spec_mat<-spec_mat_comp(weight_func,L,Lag,c_eta,lag_mat)$spec_mat
  #spec_mat[,2]  spec_math-spec_mat
  structure_func_obj<-structure_func(weight_func,spec_mat)


  Gamma_structure<-structure_func_obj$Gamma_structure
  spec_mat_structure<-structure_func_obj$spec_mat_structure
  Gamma_structure_diff<-structure_func_obj$Gamma_structure_diff
  spec_mat_structure_diff<-structure_func_obj$spec_mat_structure_diff

  # weighting of amplitude function in stopband
  omega_Gamma<-as.integer(cutoff*K/pi)
  if ((K-omega_Gamma+1)>0)
  {
    # Replicate results in McElroy-Wildi (trilemma paper)
    if (lin_eta)
    {
      eta_vec<-c(rep(1,omega_Gamma),1+rep(eta,K-omega_Gamma+1))
    } else
    {
      if (c_eta)
      {
        eta_vec<-c(rep(1,omega_Gamma+1),(1:(K-omega_Gamma))*pi/K +1)^(eta/10)
      } else
      {
        eta_vec<-c(rep(1,omega_Gamma),(1:(K-omega_Gamma+1))^(eta/2))
      }
    }
    weight_h<-weight_func*eta_vec
  } else
  {
    eta_vec<-rep(1,K+1)
    weight_h<-weight_func* eta_vec
  }
  # Frequency zero receives half weight
  weight_h[1,]<-weight_h[1,]*ifelse(c_eta,1,1/sqrt(2))
  # DFT target variable
  weight_target<-weight_h[,1]
  # Rotate all DFT's such that weight_target is real (rotation does not alter mean-square error). Note that target of strutural
  # additional structural elements and of original MDFA are the same i.e. rotation must be the same and weight_target are identical too.
  weight_h<-weight_h*exp(-1.i*Arg(weight_target))
  weight_target<-Re(weight_target*exp(-1.i*Arg(weight_target)))

  # If structure (forecasting) is imposed then we have to undo the customization weighting due to eta (no emphasis of stopband)
  if (sum(abs(weight_structure))>0)
  {
    weight_target_structure<-weight_target_structure_diff<-weight_target/eta_vec
  } else
  {
    weight_target_structure<-weight_target_structure_diff<-rep(0,K+1)
  }


  # DFT's explaining variables: target variable can be an explaining variable too
  weight_h_exp<-as.matrix(weight_h[,2:(dim(weight_h)[2])])

  # The spectral matrix is inflated in stopband: effect of eta
  spec_mat<-t(t(spec_mat)*eta_vec) #dim(spec_mat)  as.matrix(Im(spec_mat[150,]))

  # Compute design matrix and regularization matrix
  mat_obj<-mat_func(i1,i2,L,weight_h_exp,lambda_decay,lambda_cross,lambda_smooth,Lag,weight_constraint,shift_constraint,grand_mean,b0_H0,c_eta,lag_mat)

  des_mat<-mat_obj$des_mat
  reg_mat<-mat_obj$reg_mat
  reg_xtxy<-mat_obj$reg_xtxy
  w_eight<-mat_obj$w_eight

  # Solve estimation problem
  mat_x<-des_mat%*%spec_mat#spec_mat[,1]<-1
  X_new<-t(Re(mat_x))+sqrt(1+Gamma*lambda)*1.i*t(Im(mat_x))
  # xtx can be written either in Re(t(Conj(X_new))%*%X_new) or as below:
  xtx<-t(Re(X_new))%*%Re(X_new)+t(Im(X_new))%*%Im(X_new)
  # If abs(weight_structure)>0 then additional structure is instilled into MDFA-criterion: one-step ahead MSE
  if (sum(abs(weight_structure))>0)
  {
    X_new_structure<-t(abs(weight_structure[1])*des_mat%*%spec_mat_structure)
    xtx_structure<-t(Re(X_new_structure))%*%Re(X_new_structure)+t(Im(X_new_structure))%*%Im(X_new_structure)

    X_new_structure_diff<-t(abs(weight_structure[2])*des_mat%*%spec_mat_structure_diff)
    xtx_structure_diff<-t(Re(X_new_structure_diff))%*%Re(X_new_structure_diff)+t(Im(X_new_structure_diff))%*%Im(X_new_structure_diff)

    Gamma_structure_diff<-abs(weight_structure)[2]*Gamma_structure
    Gamma_structure<-abs(weight_structure)[1]*Gamma_structure

  } else
  {
    X_new_structure<-x_new_structure_diff<-0*des_mat%*%spec_mat_structure
    xtx_structure<-xtx_structure_diff<-0*t(Re(X_new_structure))%*%Re(X_new_structure)+t(Im(X_new_structure))%*%Im(X_new_structure)
  }
  # The filter restrictions (i1<-T and/or i2<-T) appear as constants on the right hand-side of the equation:
  xtxy<-t(Re(t(w_eight)%*%spec_mat)%*%Re(t(spec_mat)%*%t(des_mat))+
            Im(t(w_eight)%*%t(t(spec_mat)*sqrt(1+Gamma*lambda)))%*%Im(t(t(t(spec_mat)*sqrt(1+Gamma*lambda)))%*%t(des_mat)))
  # scaler makes scales of regularization and unconstrained optimization `similar'
  scaler<-mean(diag(xtx))




  if (sum(abs(weight_structure))>0)
  {
    X_inv<-solve(xtx+xtx_structure+xtx_structure_diff+scaler*reg_mat)
    bh<-as.vector(X_inv%*%(((t(Re(X_new)*weight_target))%*%Gamma)+
                             ((t(Re(X_new_structure)*weight_target_structure))%*%Gamma_structure)+
                             ((t(Re(X_new_structure_diff)*weight_target_structure_diff))%*%Gamma_structure_diff)-
                             xtxy-scaler*reg_xtxy))
  } else
  {
    X_inv<-solve(xtx+scaler*reg_mat)
    bh<-as.vector(X_inv%*%(((t(Re(X_new)*weight_target))%*%Gamma)-xtxy-scaler*reg_xtxy)) #weight_target[1]<-1
  }

  b<-matrix(nrow=L,ncol=length(weight_h_exp[1,]))
  # Reconstruct original parameters from the set of possibly contrained ones
  bhh<-t(des_mat)%*%bh
  for (k in 1:L)
  {
    b[k,]<-bhh[(k)+(0:(length(weight_h_exp[1,])-1))*L]
  }



  weight_cm<-matrix(w_eight,ncol=(length(weight_h_exp[1,])))
  # Add level and/or time-shift constraints (if i1<-F and i2<-F then this matrix is zero)
  b<-b+weight_cm


  # Transferfunction
  trffkt<-matrix(nrow=K+1,ncol=length(weight_h_exp[1,]))
  trffkth<-trffkt
  trffkt[1,]<-apply(b,2,sum)
  trffkth[1,]<-trffkt[1,]


  for (j in 1:length(weight_h_exp[1,]))
  {
    for (k in 0:(K))
    {
      trffkt[k+1,j]<-(b[,j]%*%exp(1.i*k*(0:(L-1))*pi/(K)))
    }
  }


  trt<-apply(((trffkt)*exp(1.i*(0-Lag)*pi*(0:(K))/K))*weight_h_exp,1,sum)
  # DFA criterion which accounts for customization but not for regularization term
  # MDFA-Legacy : new normalization for all error terms below
  rever<-sum(abs(Gamma*weight_target-Re(trt)-1.i*sqrt(1+lambda*Gamma)*Im(trt))^2)*pi/(K+1)
  # MS-filter error : DFA-criterion without effects by lambda or eta (one must divide spectrum by eta_vec)
  MS_error<-sum((abs(Gamma*weight_target-trt)/eta_vec)^2)*pi/(K+1)
  # Definition of Accuracy, time-shift and noise suppression terms
  Gamma_cp<-Gamma[1+0:as.integer(K*(cutoff/pi))]
  Gamma_cn<-Gamma[(2+as.integer(K*(cutoff/pi))):(K+1)]
  trt_cp<-(trt/eta_vec)[1+0:as.integer(K*(cutoff/pi))]
  trt_cn<-(trt/eta_vec)[(2+as.integer(K*(cutoff/pi))):(K+1)]
  weight_target_cp<-(weight_target/eta_vec)[1+0:as.integer(K*(cutoff/pi))]
  weight_target_cn<-(weight_target/eta_vec)[(2+as.integer(K*(cutoff/pi))):(K+1)]
  # define singular observations
  Accuracy<-sum(abs(Gamma_cp*weight_target_cp-abs(trt_cp))^2)*2*pi/(K+1)
  Timeliness<-4*sum(abs(Gamma_cp)*abs(trt_cp)*sin(Arg(trt_cp)/2)^2*weight_target_cp)*2*pi/(K+1)
  Smoothness<-sum(abs(Gamma_cn*weight_target_cn-abs(trt_cn))^2)*2*pi/(K+1)
  Shift_stopband<-4*sum(abs(Gamma_cn)*abs(trt_cn)*sin(Arg(trt_cn)/2)^2*weight_target_cn)*2*pi/(K+1)

# The following computations are time-consuming if K is `large'
# By default they are skipped i.e. troikaner=F
# Additional output: degrees of freedom, AIC
  if (troikaner)
  {
    # The following derivations of the DFA-criterion are equivalent
    # They are identical with rever (up to normalization by (2*(K+1)^2))
    trth<-((X_new)%*%(X_inv%*%t(Re(X_new))))%*%(weight_target*Gamma)
    # The projection matrix
    Proj_mat<-((X_new)%*%(X_inv%*%t(Re(X_new))))
    # The residual projection matrix
    res_mat<-diag(rep(1,dim(Proj_mat)[1]))-Proj_mat
    # DFA criterion: first possibility (all three variants are identical)
    sum(abs(res_mat%*%(weight_target*Gamma))^2)*pi/(K+1)
    # Residuals (DFT of filter error)
    resi<-res_mat%*%(weight_target*Gamma)
    # DFA criterion: second possibility
    t(Conj(resi))%*%resi*pi/(K+1)
    t((weight_target*Gamma))%*%(t(Conj(res_mat))%*%(res_mat))%*%(weight_target*Gamma)*pi/(K+1)
    degrees_freedom<-2*Re(sum(diag(t(Conj(res_mat))%*%(res_mat))))
    #  Re(t(Conj(res_mat))%*%(res_mat))-Re(res_mat)
    freezed_degrees<-2*K+1-degrees_freedom
    # The following expression equates to zero i.e. projection matrix is a projection for real part
    #   Re(t(Conj(res_mat))%*%(res_mat))-Re(res_mat)

    M<-((X_new)%*%(X_inv%*%t(Conj(X_new))))
    res_M<-diag(rep(1,dim(M)[1]))-M
    # The following is a real number but due to numerical rounding errors we prefer to take the real part of the result
    freezed_degrees_new<-Re(K+1-sum(eigen(res_M)$values))
    # Troikaner
    tic<-ifelse(degrees_freedom<2*K+1&degrees_freedom>1,log(rever)+2*(K-degrees_freedom+1)/(degrees_freedom-2),NA)
    tic<-log(rever)+2*freezed_degrees_new/(2*K)
    return(list(b=b,trffkt=trffkt,rever=rever,degrees_freedom=degrees_freedom,tic=tic,freezed_degrees=freezed_degrees,Accuracy=Accuracy,Smoothness=Smoothness,Timeliness=Timeliness,MS_error=MS_error,freezed_degrees_new=freezed_degrees_new))
  } else
  {
# Simplified (shorter) return
    return(list(b=b,trffkt=trffkt,rever=rever,Accuracy=Accuracy,Smoothness=Smoothness,Timeliness=Timeliness,MS_error=MS_error))
  }
}









#' This function sets-up the spectral matrix for closed-form solution of the optimization problem in the multivariate case: it relies on the multivariate DFT; it organizes lead/lags in the case of mixed-frequency data.
#'
#' @param weight_func DFT-matrix or alternative (for example model-based) estimate: first column is the target variable, additional columns are explanatory variables
#' @param L Filter-length
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param c_eta Boolean weight of frequency zero (default is F)
#' @param lag_mat Matrix for implementing effective lags in a mixed-frequency setting
#' @return spec_mat Spectral matrix set-up for closed-form solution of optimization problem
#' @export
#'

spec_mat_comp<-function(weight_func,L,Lag,c_eta,lag_mat)
{
  K<-length(weight_func[,1])-1
  weight_h<-weight_func
  # Frequency zero receives half weight
  # Chris mod: avoid halving
# MDFA-Legacy: the DFT in frequency zero is weighted by 1/sqrt(2) (the weight of frequency zero is 0.5 in mean-square)
  if (!c_eta)
    weight_h[1,]<-weight_h[1,]/sqrt(2)
  # Extract DFT target variable (first column)
  weight_target<-weight_h[,1]
  # Rotate all DFT's such that weight_target is real (rotation does not alter mean-square error)
  weight_h<-weight_h*exp(-1.i*Arg(weight_target))#Im(exp(-1.i*Arg(weight_target)))
  weight_target<-weight_target*exp(-1.i*Arg(weight_target))
  # DFT's explaining variables (target variable can be an explaining variable too)
  weight_h_exp<-as.matrix(weight_h[,2:(dim(weight_h)[2])])
  spec_mat<-as.vector(t(as.matrix(weight_h_exp[1,])%*%t(as.matrix(rep(1,L)))))

  for (j in 1:(K))#j<-1  h<-2  lag_mat<-matrix(rep(0:1,3),ncol=3)   Lag<-2
  {
    omegak<-j*pi/K
# We feed the daily structure of the mixed_frequency data lag_mat
    exp_mat<-matrix(nrow=L,ncol=ncol(weight_h_exp))
    for (h in 1:ncol(exp_mat))
      exp_mat[,h]<-exp(1.i*omegak*(lag_mat[,h+1]-Lag))
# Inclusion of the time-varying lag-structure of the mixed-frequency data
    spec_mat<-cbind(spec_mat,as.vector(t(weight_h_exp[j+1,]*t(exp_mat))))
  }
  dim(spec_mat)
  return(list(spec_mat=spec_mat))#as.matrix(Re(spec_mat[1,]))
}











#' This function sets-up the design matrix of the generic optimization problem: it combines data, regularization features and filter constraints
#'
#' @param i1 Boolean. If T a first-order filter constraint in frequency zero is obtained: amplitude of real-time filter must match weight_constraint (handles integration order one)
#' @param i2 Boolean. If T a second-order filter constraint in frequency zero is obtained: time-shift of real-time filter must match target (together with i1 handles integration order two)
#' @param L Filter-length
#' @param weight_h_exp DFT of explanatory variables
#' @param lambda_decay Regularization: decay term
#' @param lambda_cross Regularization: cross-sectional term
#' @param lambda_smooth Regularization: smoothness term
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param weight_constraint Constraint vector in the case i1==T
#' @param shift_constraint Constraint vector in the case i2==T
#' @param grand_mean Boolean: if T then a grand-mean parametrization is imposed (default is F)
#' @param b0_H0 Regularization: shrinkage target (arbitrary designs can be replicated by imposing strong regularization)
#' @param c_eta Boolean: impose mild/strong smoothness customization (default is F)
#' @param lag_mat Matrix for implementing effective lags in a mixed-frequency setting
#' @return des_mat Design matrix
#' @return reg_mat Bilinear form (regularization matrix)
#' @return reg_xtxy Constants from regularization
#' @return w_eight Constants from filter constraints
#'

mat_func<-function(i1,i2,L,weight_h_exp,lambda_decay,lambda_cross,lambda_smooth,Lag,weight_constraint,shift_constraint,grand_mean,b0_H0,c_eta,lag_mat)
{



# MDFA-Legacy: new functions
# Regularization
  reg_mat_obj<-reg_mat_func(weight_h_exp,L,c_eta,lambda_decay,Lag,grand_mean,lambda_cross,lambda_smooth,lag_mat)

  Q_smooth<-reg_mat_obj$Q_smooth
  Q_cross<-reg_mat_obj$Q_cross
  Q_decay<-reg_mat_obj$Q_decay
  lambda_decay<-reg_mat_obj$lambda_decay
  lambda_smooth<-reg_mat_obj$lambda_smooth
  lambda_cross<-reg_mat_obj$lambda_cross

# MDFA-Legacy: new functions
# weight vector for constraints
  w_eight<-w_eight_func(i1,i2,Lag,weight_constraint,shift_constraint,L,weight_h_exp)$w_eight

# MDFA-Legacy: new functions
# Grand-mean and Constraints

  des_mat<-des_mat_func(i2,i1,L,weight_h_exp,weight_constraint,shift_constraint,Lag)$des_mat

# Transforming back from grand-mean if necessary (des_mat assumes grand-mean)
  if (!grand_mean&length(weight_h_exp[1,])>1)
  {
    Q_centraldev_original<-centraldev_original_func(L,weight_h_exp)$Q_centraldev_original
  } else
  {
    Q_centraldev_original<-NULL
  }
  if (!grand_mean&length(weight_h_exp[1,])>1)
  {
# apply A^{-1} to A*R giving R where A and R are defined in MDFA_Legacy book
#and des_mat=A*R i.e.  t(Q_centraldev_original%*%t(des_mat)) is R (called des_mat in my code)
    des_mat<-t(Q_centraldev_original%*%t(des_mat))
  }

# Here we fold all three regularizations (cross, smooth and decay) into a single reg-matrix
  if ((length(weight_h_exp[1,])>1))
  {
    if (grand_mean)
    {
      reg_t<-(Q_smooth+Q_decay+t(Q_centraldev_original)%*%Q_cross%*%Q_centraldev_original)
    } else
    {
      reg_t<-(Q_smooth+Q_decay+Q_cross)
    }
  } else
  {
    reg_t<-(Q_smooth+Q_decay)
  }



# Normalize regularization terms (which are multiplied by des_mat) in order to disentangle i1/i2 effects
  reg_mat<-(des_mat)%*%reg_t%*%t(des_mat)#sum(diag(reg_mat))
  if (is.null(b0_H0))
    b0_H0<-rep(0,length(w_eight))
  b0_H0<-as.vector(b0_H0)
  if (lambda_smooth+lambda_decay[2]+lambda_cross>0)
  {
    disentangle_des_mat_effect<-sum(diag(reg_t))/sum(diag(reg_mat))
    if (abs(sum(diag(reg_mat)))>0)
    {
      disentangle_des_mat_effect<-sum(diag(reg_t))/sum(diag(reg_mat))
    } else
    {
      disentangle_des_mat_effect<-1
    }

    reg_mat<-reg_mat*disentangle_des_mat_effect
    reg_xtxy<-des_mat%*%reg_t%*%(w_eight-b0_H0)*(disentangle_des_mat_effect)#+t(w_eight)%*%reg_t%*%t(des_mat)
  } else
  {
    reg_xtxy<-des_mat%*%reg_t%*%(w_eight-b0_H0)#+t(w_eight)%*%reg_t%*%t(des_mat)
    dim(des_mat)
    dim(reg_t)
  }
  return(list(des_mat=des_mat,reg_mat=reg_mat,reg_xtxy=reg_xtxy,w_eight=w_eight))
}









#' This function transforms grand-mean parametrization back to original coefficients
#'
#' @param L Filter-length
#' @param weight_h_exp DFT of explanatory variables
#' @return Q_centraldev_original Transformation matrix
#'

centraldev_original_func<-function(L,weight_h_exp)
{
# Grand-mean parametrization: des_mat is implemented in terms of grand-mean and Q_centraldev_original can
# be used to invert back to original coefficients
  Q_centraldev_original<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))


  Q_centraldev_original<-diag(rep(1,L*length(weight_h_exp[1,])))
  if (L>1)
  {
    diag(Q_centraldev_original[1:L,L+1:L])<-rep(-1,L)
    for (i in 2:length(weight_h_exp[1,]))   #i<-2
    {
      diag(Q_centraldev_original[(i-1)*L+1:L,1:L])<-rep(1,L)
      diag(Q_centraldev_original[(i-1)*L+1:L,(i-1)*L+1:L])<-rep(1,L)
      diag(Q_centraldev_original[1:L,(i-1)*L+1:L])<-rep(-1,L)
    }
  }
  Q_centraldev_original<-solve(Q_centraldev_original)
  return(list(Q_centraldev_original=Q_centraldev_original))
}








#' This function calls the regularization settings (bilinear forms) for cross-correlation, smoothness and decay patterns
#'
#' @param weight_h_exp DFT of explanatory variables
#' @param L Filter-length
#' @param c_eta Boolean: impose mild/strong smoothness customization (default is F)
#' @param lambda_decay Regularization: decay term
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param grand_mean Boolean: if T then a grand-mean parametrization is imposed (default is F)
#' @param lambda_cross Regularization: cross-sectional term
#' @param lambda_smooth Regularization: smoothness term
#' @param lag_mat Matrix for implementing effective lags in a mixed-frequency setting
#' @return Q_smooth Smoothness term (matrix)
#' @return Q_cross Cross-sectional term (matrix)
#' @return Q_decay Decay term (matrix)
#' @return lambda_decay Re-parameterized decay-weight
#' @return lambda_smooth Re-parameterized smoothness-weight
#' @return lambda_cross Re-parameterized cross-sectional-weight
#'

reg_mat_func<-function(weight_h_exp,L,c_eta,lambda_decay,Lag,grand_mean,lambda_cross,lambda_smooth,lag_mat)
{
  lambda_smooth<-100*tan(min(abs(lambda_smooth),0.999999999)*pi/2)
  lambda_cross<-100*tan(min(abs(lambda_cross),0.9999999999)*pi/2)
  if (c_eta)
  {
    # Chris replication: imposing non-linear transform to lambda_decay[1] too
    lambda_decay<-c(min(abs(tan(min(abs(lambda_decay[1]),0.999999999)*pi/2)),1),100*tan(min(abs(lambda_decay[2]),0.999999999)*pi/2))
  } else
  {
    lambda_decay<-c(min(abs(lambda_decay[1]),1),100*tan(min(abs(lambda_decay[2]),0.999999999)*pi/2))
  }
  # New regularization function

  Q_obj<-Q_reg_func(L,weight_h_exp,lambda_decay,lambda_smooth,lambda_cross,Lag,lag_mat,grand_mean,c_eta)

  Q_smooth<-Q_obj$Q_smooth
  Q_decay<-Q_obj$Q_decay
  Q_cross<-Q_obj$Q_cross
  # The Q_smooth and Q_decay matrices address regularizations for original unconstrained parameters (Therefore dimension L^2)
  # At the end, the matrix des_mat is used to map these regularizations to central-deviance parameters
  # accounting for first order constraints!

  return(list(Q_smooth=Q_smooth,Q_cross=Q_cross,Q_decay=Q_decay,lambda_decay=lambda_decay,
              lambda_smooth=lambda_smooth,lambda_cross=lambda_cross))
}






#' This function sets-up the regularization matrices
#'
#' @param L Filter-length
#' @param weight_h_exp DFT of explanatory variables
#' @param lambda_decay Regularization: decay term
#' @param lambda_smooth Regularization: smoothness term
#' @param lambda_cross Regularization: cross-sectional term
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param lag_mat Matrix for implementing effective lags in a mixed-frequency setting
#' @param grand_mean Boolean: if T then a grand-mean parametrization is imposed (default is F)
#' @param c_eta Boolean: impose mild/strong smoothness customization (default is F)
#' @return Q_smooth Smoothness term (matrix)
#' @return Q_decay Decay term (matrix)
#' @return Q_cross Cross-sectional term (matrix)
#'

# Function for setting-up the bilinear forms of the regularization
# It is a technical function which does not deserve specific description

Q_reg_func<-function(L,weight_h_exp,lambda_decay,lambda_smooth,lambda_cross,Lag,lag_mat,grand_mean,c_eta)
{
  Q_smooth<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  Q_decay<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  # Cross-sectional regularization if dimension>1
  if ((length(weight_h_exp[1,])>1))
  {
    # The cross-sectional regularization is conveniently implemented on central-deviance parameters. The regularization is expressed on the
    # unconstrained central-deviance parameters (dimension L), then mapped to the original (unconstrained) parameters (dimension L) with Q_centraldev_original
    # and then maped back to central-deviance with constraint (dim L-1) with des_mat (mathematically unnecessarily complicate but more convenient to implement in code).
    Q_cross<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
    Q_centraldev_original<-matrix(data=rep(0,((L)*length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
  } else
  {
    # 16.08.2012
    Q_cross<-NULL
  }
  for (i in 1:L)
  {
    # For symmetric filters or any historical filter with Lag>0 the decay must be symmetric about b_max(0,Lag)
    # lambda_decay is a 2-dim vector: the first component controls for the exponential decay and the second accounts for the strength of the regularization
    # The maximal weight is limited to 1e+4
    Q_decay[i,i]<-min((1+lambda_decay[1])^(2*abs(i-1-max(0,Lag))),1e+4)
    # Original idea: with mixed-frequency data the decay should account for the effective lag (of the lower frequency data) as measured on the high-frequency scale: this is not a clever idea...
    #   Q_decay[i,i]<-min((1+lambda_decay[1])^(2*abs(lag_mat[i,1+1]-max(0,Lag))),1e+4)

    if (L>4)
    {
      if(i==1)
      {
        Q_smooth[i,i:(i+2)]<-c(1,-2,1)
      } else
      {
        if(i==2)
        {
          Q_smooth[i,(i-1):(i+2)]<-c(-2,5,-4,1)
        } else
        {
          if(i==L)
          {
            Q_smooth[i,(i-2):i]<-c(1,-2,1)
          } else
          {
            if(i==L-1)
            {
              Q_smooth[i,(i-2):(i+1)]<-c(1,-4,5,-2)
            } else
            {
              Q_smooth[i,(i-2):(i+2)]<-c(1,-4,6,-4,1)
            }
          }
        }
      }
    } else
    {
      print("L<=4: no smoothness regularization imposed!!!!!!!!!")
    }
  }

  if (length(weight_h_exp[1,])>1)
  {

    for (j in 1:max(1,(length(weight_h_exp[1,])-1)))   #j<-1
    {
      if (L>4)
        Q_smooth[j*L+1:L,j*L+1:L]<-Q_smooth[1:L,1:L]
      for (i in 1:L)
        Q_decay[j*L+i,j*L+i]<-min((1+lambda_decay[1])^(2*abs(i-1-max(0,Lag))),1e+4)
# Original idea: with mixed-frequency data the decay should account for the effective lag (of the lower frequency data) as measured on the high-frequency scale: this is not a clever idea...
#      Q_decay[j*L+i,j*L+i]<-min((1+lambda_decay[1])^(2*abs(lag_mat[i,1+j+1]-max(0,Lag))),1e+4)

    }
    Q_centraldev_original<-diag(rep(1,L*length(weight_h_exp[1,])))
    if (L>1)
    {
      diag(Q_centraldev_original[1:L,L+1:L])<-rep(-1,L)
      for (i in 2:length(weight_h_exp[1,]))   #i<-2
      {
        diag(Q_centraldev_original[(i-1)*L+1:L,1:L])<-rep(1,L)
        diag(Q_centraldev_original[(i-1)*L+1:L,(i-1)*L+1:L])<-rep(1,L)
        diag(Q_centraldev_original[1:L,(i-1)*L+1:L])<-rep(-1,L)
      }
    }
    Q_centraldev_original<-solve(Q_centraldev_original)

    # 06.08.2012: the following if allows for either grand-mean parametrization (I-MDFA version prior to 30.07.2012) or original parameters
    #   If grand_mean==T then the code replicates I-MDFA as released prior 30.07.2012
    #   If grand_mean==F then the new parametrization is used.
    #   Differences between both approaches: see section 7.2 of my elements paper posted on SEFBlog (both codes are identical when no regularization is imposed. Otherwise the later version (grand_mean==F) is logically more consistent becuase it treats all series identically (no asymmetry)).
    if (grand_mean)
    {
      diag(Q_cross[L+1:((length(weight_h_exp[1,])-1)*L),L+1:((length(weight_h_exp[1,])-1)*L)])<-
        rep(1,((length(weight_h_exp[1,])-1)*L))
    } else
    {
      #30.07.2012:new definition (parametrization) of Q_cross (Lambda_{cross} in the elements-paper)
      diag(Q_cross)<-1
      for (i in 1:length(weight_h_exp[1,]))
      {
        for (j in 1:L)
        {
          Q_cross[(i-1)*L+j,j+(0:(length(weight_h_exp[1,])-1))*L]<-Q_cross[(i-1)*L+j,j+(0:(length(weight_h_exp[1,])-1))*L]-1/length(weight_h_exp[1,])
        }
      }
    }
  } else
  {
    # define matrix for univariate case
    Q_centraldev_original<-NULL
  }
  # Normalizing the troika are new: disentangle the effect by L
  Q_decay<-Q_decay*lambda_decay[2]
  Q_cross<-Q_cross*lambda_cross                   #Qh<-Q_cross
  Q_smooth<-Q_smooth*lambda_smooth
  if (lambda_decay[2]>0)
  {
    # The second parameter in lambda_decay accounts for the strength of the regularization
    Q_decay<-lambda_decay[2]*(Q_decay/(sum(diag(Q_decay))))
  }
  if (lambda_cross>0)
  {
    if (c_eta)
    {
      Q_cross<-lambda_cross^2*(Q_cross/(sum(diag(Q_cross))))
    } else
    {
      Q_cross<-lambda_cross*(Q_cross/(sum(diag(Q_cross))))
    }
  }
  if (lambda_smooth>0&L>4)
  {
    Q_smooth<-lambda_smooth*(Q_smooth/(sum(diag(Q_smooth))))
  }
  return(list(Q_smooth=Q_smooth,Q_decay=Q_decay,Q_cross=Q_cross))
}







#' This function specifies the constant in case of filter constraints with/without regularization
#'
#' @param i1 Boolean. If T a first-order filter constraint in frequency zero is obtained: amplitude of real-time filter must match weight_constraint (handles integration order one)
#' @param i2 Boolean. If T a second-order filter constraint in frequency zero is obtained: time-shift of real-time filter must match target (together with i1 handles integration order two)
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @param weight_constraint Constraint vector in the case i1==T
#' @param shift_constraint Constraint vector in the case i2==T
#' @param L Filter-length
#' @param weight_h_exp DFT of explanatory variables
#' @return w_eight Vector of filter weights entering closed-form solution
#'

w_eight_func<-function(i1,i2,Lag,weight_constraint,shift_constraint,L,weight_h_exp)
{
  if (i1)
  {
    if (i2)
    {
#               impose constraints to b_Lag, b_{Lag+1} instead of b_{L-1} and b_L
#               Therefore the decay regularization does not potentially conflict with filter constraints
      if (Lag<1)
      {
# corrected time-shift expression if A(0) different from 1
        w_eight<-c(-(Lag-1)*weight_constraint[1]-shift_constraint[1]*weight_constraint[1],
                   Lag*weight_constraint[1]+shift_constraint[1]*weight_constraint[1],rep(0,L-2))
      } else
      {
# corrected time-shift expression if A(0) different from 1
        w_eight<-c(rep(0,Lag),weight_constraint[1]-shift_constraint[1]*weight_constraint[1],
                   shift_constraint[1]*weight_constraint[1],rep(0,L-Lag-2))
      }

      if (length(weight_h_exp[1,])>1)
      {
        for (j in 2:length(weight_h_exp[1,]))
        {

          if (Lag<1)
          {
# corrected time-shift expression if A(0) different from 1
            w_eight<-c(w_eight,-(Lag-1)*weight_constraint[j]-shift_constraint[j]*weight_constraint[j],
                       Lag*weight_constraint[j]+shift_constraint[j]*weight_constraint[j],rep(0,L-2))
          } else
          {
# corrected time-shift expression if A(0) different from 1
            w_eight<-c(w_eight,c(rep(0,Lag),weight_constraint[j]-shift_constraint[j]*weight_constraint[j],
                                 shift_constraint[j]*weight_constraint[j],rep(0,L-Lag-2)))
          }
        }
      }
    } else
    {
      if (Lag<1)
      {
        w_eight<-c(weight_constraint[1],rep(0,L-1))
      } else
      {
        w_eight<-c(rep(0,Lag),weight_constraint[1],rep(0,L-Lag-1))
      }
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 2:length(weight_h_exp[1,]))
        {
          if (Lag<1)
          {
            w_eight<-c(w_eight,weight_constraint[j],rep(0,L-1))
          } else
          {
            w_eight<-c(w_eight,rep(0,Lag),weight_constraint[j],rep(0,L-Lag-1))
          }
        }
      }
    }
  } else
  {
# MDFA-Legacy : the case i2==T i1==F is fixed, at last.
    if (i2)
    {
      if (Lag<1)
      {
        w_eight<-rep(0,L)
      } else
      {
        w_eight<-rep(0,L)
      }

      if (length(weight_h_exp[1,])>1)
      {
        for (j in 2:length(weight_h_exp[1,]))
        {
          if (Lag<1)
          {
            w_eight<-c(w_eight,rep(0,L))
          } else
          {
            w_eight<-c(w_eight,rep(0,L))
          }
        }
      }
    } else
    {
      w_eight<-rep(0,L*length(weight_h_exp[1,]))
    }
  }
  return(list(w_eight= w_eight))

}





#' This function sets-up the data-based part of the design matrix, accounting for i1/i2 constraints (this is the effective design-matrix in the absence of regularization). It is originally written in grand-mean parametrization (which makes it more tedious).  Singularities in the case i1=F, i2=T are avoided by `displacing' the constraints by a small number.
#'
#' @param i2 Boolean. If T a second-order filter constraint in frequency zero is obtained: time-shift of real-time filter must match target (together with i1 handles integration order two)
#' @param i1 Boolean. If T a first-order filter constraint in frequency zero is obtained: amplitude of real-time filter must match weight_constraint (handles integration order one)
#' @param L Filter-length
#' @param weight_h_exp DFT of explanatory variables
#' @param weight_constraint Constraint vector in the case i1==T
#' @param shift_constraint Constraint vector in the case i2==T
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @return des_mat Data-based design matrix
#'

des_mat_func<-function(i2,i1,L,weight_h_exp,weight_constraint,shift_constraint,Lag)
{
  # Here we implement the matrix which links freely determined central-deviance parameters and constrained original parameters
  # In the MDFA_Legacy book t(des_mat) corresponds to A%*%R (R links constrained and unconstrained parameters and A maps central-deviance to original parameters)
  #   Please note that:
  #   1. Here I'm working with central-deviance parameters
  #   2. The same matrix R applies to either parameter set
  #   3. If I work with central-deviance parameters then R maps the freely determined set to the constrained (central-deviance)
  #       and A then maps the constrained (central-deviance) set to original constrained parameters.

  if (i2)
  {
    if (i1)
    {
      # First and second order restrictions
      des_mat<-matrix(data=rep(0,(L-2)*L*(length(weight_h_exp[1,]))^2),nrow=(L-2)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))

      for (i in 1:(L-2))
      {
        if (Lag<1)
        {
          des_mat[i,i+2+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,1+(0:(length(weight_h_exp[1,])-1))*L]<-i
          des_mat[i,2+(0:(length(weight_h_exp[1,])-1))*L]<--(i+1)

        } else
        {
          des_mat[i,ifelse(i<Lag+1,i,i+2)+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,Lag+1+(0:(length(weight_h_exp[1,])-1))*L]<-ifelse(i<Lag+1,-(Lag+2-i),i-Lag)
          des_mat[i,Lag+2+(0:(length(weight_h_exp[1,])-1))*L]<-ifelse(i<Lag+1,(Lag+1-i),-(i-Lag+1))
        }
      }
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 1:max(1,(length(weight_h_exp[1,])-1)))
        {
          for (i in 1:(L-2))
          {

            if (Lag<1)
            {
              des_mat[i+j*(L-2),i+2]<--1
              des_mat[i+j*(L-2),1]<--i
              des_mat[i+j*(L-2),2]<-(i+1)
            } else
            {
              des_mat[i+j*(L-2),ifelse(i<Lag+1,i,i+2)]<--1
              des_mat[i+j*(L-2),Lag+1]<--ifelse(i<Lag+1,-(Lag+2-i),i-Lag)
              des_mat[i+j*(L-2),Lag+2]<--ifelse(i<Lag+1,(Lag+1-i),-(i-Lag+1))
            }

            if (Lag<1)
            {
              des_mat[i+j*(L-2),i+2+j*L]<-1
              des_mat[i+j*(L-2),1+j*L]<-i
              des_mat[i+j*(L-2),2+j*L]<--(i+1)
            } else
            {
              des_mat[i+j*(L-2),ifelse(i<Lag+1,i+j*L,i+2+j*L)]<-1
              des_mat[i+j*(L-2),Lag+1+j*L]<-ifelse(i<Lag+1,-(Lag+2-i),i-Lag)
              des_mat[i+j*(L-2),Lag+2+j*L]<-ifelse(i<Lag+1,(Lag+1-i),-(i-Lag+1))
            }
          }
        }
      }
    } else
    {

      des_mat<-matrix(data=rep(0,(L-1)*L*(length(weight_h_exp[1,]))^2),nrow=(L-1)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
# MDFA-Legacy : the case i2==T i1==F is fixed, at last
# Avoid singularities in the time-shift specification
      epsi<-1.e-02
      if (Lag>=1)
      shift_constraint[which(abs(shift_constraint)<epsi)]<-shift_constraint[which(abs(shift_constraint)<epsi)]+epsi
      if (Lag<1)
        shift_constraint[abs(-Lag-shift_constraint)<epsi]<-shift_constraint[which(abs(-Lag-shift_constraint)<epsi)]+epsi
      for (i in 1:(L-1))
      {
        if (Lag<1)
        {
          # Forecast and nowcast
          # we have to differentiate the case with vanishing shift and non-vansihing shift
          #    For a vanishing shift we select b2 for the constraint because then we are sure that the equations can be solved
          #    For a non-vanishing shift we select b1 because then we are sure that the equations can be solved (irrespective of the shift)
          # See MDFA-legacy book for reference

            # shift is different from zero: b1 is isolated (new formulas)
          des_mat[i,ifelse(i<1,i,i+1)+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,1+(0:(length(weight_h_exp[1,])-1))*L]<--(-Lag+i-shift_constraint[1])/
              (-Lag-shift_constraint[1])
        } else
        {
          # Backcast: we have to differentiate the case with vanishing shift and non-vansihing shift
          #    For a vanishing shift we select b2 for the constraint because then we are sure that the equations can be solved
          #    For a non-vanishing shift we select b1 because then we are sure that the equations can be solved (irrespective of the shift)
          # See MDFA-legacy book for reference
          # shift is different from zero: b1 is isolated (new formulas: the signs change because of multiplication with -s)
          des_mat[i,ifelse(i<Lag+1,i,i+1)+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,Lag+1+(0:(length(weight_h_exp[1,])-1))*L]<-
              ifelse(i<Lag+1,(-(Lag+1-i)-shift_constraint[1])/shift_constraint[1],(-(Lag-i)-shift_constraint[1])/shift_constraint[1])

        }
      }
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 1:max(1,(length(weight_h_exp[1,])-1)))#j<-1
        {
          for (i in 1:(L-1))
          {
            if (Lag<1)
            {
              # shift is different from zero: b1 is isolated (new formulas)
              des_mat[i+j*(L-1),ifelse(i<1,i,i+1)]<--1
              des_mat[i+j*(L-1),1]<-(-Lag+i-shift_constraint[j+1])/(-Lag-shift_constraint[j+1])

            } else
            {
              # shift is different from zero: b1 is isolated (new formulas: the signs change because of multiplication with -s)
              des_mat[i+j*(L-1),ifelse(i<Lag+1,i,i+1)]<--1
              des_mat[i+j*(L-1),Lag+1]<-
                  -ifelse(i<Lag+1,(-(Lag+1-i)-shift_constraint[j+1])/shift_constraint[j+1],(-(Lag-i)-shift_constraint[j+1])/shift_constraint[j+1])

            }

            if (Lag<1)
            {
              # shift is different from zero: b1 is isolated (new formulas)
              des_mat[i+j*(L-1),ifelse(i<1,i,i+1)+j*L]<-1
              des_mat[i+j*(L-1),1+j*L]<--(-Lag+i-shift_constraint[j+1])/(-Lag-shift_constraint[j+1])

            } else
            {
              # shift is different from zero: b1 is isolated (new formulas: the signs change because of multiplication with -s)
              des_mat[i+j*(L-1),ifelse(i<Lag+1,i,i+1)+j*L]<-1
              des_mat[i+j*(L-1),Lag+1+j*L]<-
                  ifelse(i<Lag+1,(-(Lag+1-i)-shift_constraint[j+1])/shift_constraint[j+1],(-(Lag-i)-shift_constraint[j+1])/shift_constraint[j+1])

            }

          }
        }
      }

    }
  } else
  {
    if (i1)
    {
      des_mat<-matrix(data=rep(0,(L-1)*L*(length(weight_h_exp[1,]))^2),nrow=(L-1)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
      for (i in 1:(L-1))
      {
        # The i1-constraint is imposed on b_max(0,Lag) (instead of b_L ) in order to avoid a conflict with the exponential decay requirement
        if (Lag<1)
        {
          des_mat[i,i+1+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,1+(0:(length(weight_h_exp[1,])-1))*L]<--1
        } else
        {
          # Lag cannot be larger than (L-1)/2 (symmetric filter)
          des_mat[i,ifelse(i<Lag+1,i,i+1)+(0:(length(weight_h_exp[1,])-1))*L]<-1
          des_mat[i,Lag+1+(0:(length(weight_h_exp[1,])-1))*L]<--1
        }
      }
      if (length(weight_h_exp[1,])>1)
      {
        for (j in 1:max(1,(length(weight_h_exp[1,])-1)))   #j<-1
        {
          for (i in 1:(L-1))
          {
            # The i1-constraint is imposed on b_max(0,Lag) (instead of b_L ) in order to avoid a conflict with the exponential decay requirement
            if (Lag<1)
            {
              des_mat[i+j*(L-1),i+1]<--1
              des_mat[i+j*(L-1),1]<-1
              des_mat[i+j*(L-1),i+1+j*L]<-1
              des_mat[i+j*(L-1),1+j*L]<--1
            } else
            {
              # Lag cannot be larger than (L-1)/2 (symmetric filter)
              des_mat[i+j*(L-1),ifelse(i<Lag+1,i,i+1)]<--1
              des_mat[i+j*(L-1),Lag+1]<-1
              des_mat[i+j*(L-1),ifelse(i<Lag+1,i,i+1)+j*L]<-1
              des_mat[i+j*(L-1),Lag+1+j*L]<--1
            }
          }
        }
      }
    } else
    {
      des_mat<-matrix(data=rep(0,(L)*L*(length(weight_h_exp[1,]))^2),nrow=(L)*length(weight_h_exp[1,]),ncol=(L)*length(weight_h_exp[1,]))
      for (i in 1:(L))
      {
        des_mat[i,i+(0:(length(weight_h_exp[1,])-1))*L]<-1
      }

      if (length(weight_h_exp[1,])>1)
      {

        for (j in 1:max(1,(length(weight_h_exp[1,])-1)))
        {
          for (i in 1:(L))
          {
            des_mat[i+(j)*(L),i]<--1
            des_mat[i+(j)*(L),i+j*L]<-1
          }
        }
      }
    }
  }
  return(list(des_mat=des_mat))

}








#' This function decomposes the MSE-criterion into Accuracy, Smoothness and Timeliness components
#'
#' @param Gamma Generic target specification: typically symmetric Lowpass (trend) or Bandpass (cycle) filters. Highpass and anticipative allpass (forecast) can be specified too
#' @param trffkt Complex transfer function of optimal multivariate filter
#' @param weight_func DFT-matrix
#' @param cutoff Specifies start-frequency in stopband from which Smoothness is emphasized (corresponds typically to the cutoff of the lowpass target). Is not used if eta=0.
#' @param Lag Nowcast (Lag=0), Forecast (Lag<0), Backcast (Lag>0)
#' @return Accuracy Accuracy term in decomposition of MSE
#' @return Smoothness Smoothness term in decomposition of MSE
#' @return Timeliness Timeliness term in decomposition of MSE
#' @return MS_error Sample estimate of MSE: consistent estimate in the cases lambda>0 and/or eta>0
#'

MS_decomp_total<-function(Gamma,trffkt,weight_func,cutoff,Lag)
{
  if (!(length(trffkt[,1])==length(weight_func[,1])))
  {
    len_w<-min(length(trffkt[,1]),length(weight_func[,1]))
    if (length(trffkt[,1])<length(weight_func[,1]))
    {
      len_r<-(length(weight_func[,1])-1)/(length(trffkt[,1])-1)
      weight_funch<-weight_func[c(1,(1:(len_w-1))*len_r),]
      trffkth<-trffkt
    } else
    {
      len_r<-1/((length(weight_func[,1])-1)/(length(trffkt[,1])-1))
      trffkth<-trffkt[c(1,(1:(len_w-1))*len_r),]
      weight_funch<-weight_func
    }
  } else
  {
    len_w<-length(trffkt[,1])
    weight_funch<-weight_func
    trffkth<-trffkt
    Gammah<-Gamma
  }
  if (length(Gamma)>len_w)
  {
    len_r<-(length(Gamma)-1)/(len_w-1)
    Gammah<-Gamma[c(1,(1:(len_w-1))*len_r)]
  }



  weight_h<-weight_funch
  K<-length(weight_funch[,1])-1
  weight_target<-weight_h[,1]
  # Rotate all DFT's such that weight_target is real (rotation does not alter mean-square error)
  weight_h<-weight_h*exp(-1.i*Arg(weight_target))
  weight_target<-Re(weight_target*exp(-1.i*Arg(weight_target)))
  # DFT's explaining variables: target variable can be an explaining variable too
  weight_h_exp<-as.matrix(weight_h[,2:(dim(weight_h)[2])])

  trt<-apply(((trffkth)*exp(1.i*(0-Lag)*pi*(0:(K))/K))*weight_h_exp,1,sum)
  # MS-filter error : DFA-criterion without effects by lambda or eta (one must divide spectrum by eta_vec)
  MS_error<-sum((abs(Gammah*weight_target-trt))^2)/(2*(K+1)^2)
  Gamma_cp<-Gammah[1+0:as.integer(K*(cutoff/pi))]
  Gamma_cn<-Gammah[(2+as.integer(K*(cutoff/pi))):(K+1)]
  trt_cp<-trt[1+0:as.integer(K*(cutoff/pi))]
  trt_cn<-trt[(2+as.integer(K*(cutoff/pi))):(K+1)]
  weight_target_cp<-weight_target[1+0:as.integer(K*(cutoff/pi))]
  weight_target_cn<-weight_target[(2+as.integer(K*(cutoff/pi))):(K+1)]
  # define singular observations
  Accuracy<-sum(abs(Gamma_cp*weight_target_cp-abs(trt_cp))^2)/(2*(K+1)^2)
  Timeliness<-4*sum(abs(Gamma_cp)*abs(trt_cp)*sin(Arg(trt_cp)/2)^2*weight_target_cp)/(2*(K+1)^2)
  Smoothness<-sum(abs(Gamma_cn*weight_target_cn-abs(trt_cn))^2)/(2*(K+1)^2)
  Shift_stopband<-4*sum(abs(Gamma_cn)*abs(trt_cn)*sin(Arg(trt_cn)/2)^2*weight_target_cn)/(2*(K+1)^2)

  return(list(Accuracy=Accuracy,Smoothness=Smoothness,Timeliness=Timeliness,MS_error=MS_error))
}






#' This function allows for additional structure in the optimization criteria (features are not implemented yet)
#'
#' @param weight_func DFT-matrix
#' @param spec_mat Spectral matrix
#' @return Gamma_structure Target structure
#' @return spec_mat_structure Spectral matrix structure
#' @return Gamma_structure_diff Differenced target structure
#' @return spec_mat_structure Differenced spectral matrix structure
#'



structure_func<-function(weight_func,spec_mat)
{
  K<-nrow(weight_func)-1

    # We have to initialize spect_mat_structure and spec_mat_structure_diff by an arbitrary matrix of right dimensions
  spec_mat_structure<-spec_mat_structure_diff<-spec_mat
  Gamma_structure<-Gamma_structure_diff<-rep(0,K+1)

  return(list(Gamma_structure=Gamma_structure,spec_mat_structure=spec_mat_structure,
              Gamma_structure_diff=Gamma_structure_diff,spec_mat_structure_diff=spec_mat_structure_diff        ))
}



#' This function allows for particular parametrizations of DFT (not implemented yet)
#'
#' @param weight_func DFT-matrix
#' @param white_noise Boolean for flat spectrum
#' @param synchronicity Boolean for synchronous explanatory variables
#' @return weight_func Parameterized DFT
#'

white_noise_synchronicity<-function(weight_func,white_noise,synchronicity)
{

  return(list(weight_func=weight_func))
}
