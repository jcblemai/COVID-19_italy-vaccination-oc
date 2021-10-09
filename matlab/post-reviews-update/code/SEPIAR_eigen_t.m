function [R_t,e_t]=SEPIAR_eigen_t(...
            betaP,betaI,betaA,epsilon,chiE,chiP,chiI,chiA,...
            mu,deltaE,deltaP,sigma,gammaI,eta,alpha,gammaA,...
            St,Et,Pt,It,At,Rt,MS,ME,MP,MI,MA,MR,n,W,Wpi)
    % some useful matrices
    u=ones(1,n); % all-one row vector of size n
    u1n=1:n; % 1-to-n row vector
    U=sparse(u1n,u1n,u,n,n); % identity matrix of size n
    Z=sparse(zeros(n)); % all-zero n-by-n matrix
    
    % aggregated transition rates
    phiEq=mu+deltaE; phiE_v=phiEq+chiE;
    phiPq=mu+deltaP; phiP_v=phiPq+chiP;
    phiIh=mu+alpha+gammaI+eta; phiI_v=phiIh+chiI; 
    phiAq=mu+gammaA; phiA_v=phiAq+chiA;
    
    % sparse-diagonalizing distributed parameters 
    St=sparse(u1n,u1n,St,n,n); Et=sparse(u1n,u1n,Et,n,n);
    Pt=sparse(u1n,u1n,Pt,n,n); It=sparse(u1n,u1n,It,n,n);
    At=sparse(u1n,u1n,At,n,n); Rt=sparse(u1n,u1n,Rt,n,n);
    epsilon=sparse(u1n,u1n,epsilon,n,n);
    phiE=sparse(u1n,u1n,phiE_v,n,n); chiE=sparse(u1n,u1n,chiE,n,n); 
    phiP=sparse(u1n,u1n,phiP_v,n,n); chiP=sparse(u1n,u1n,chiP,n,n); 
    phiI=sparse(u1n,u1n,phiI_v,n,n); chiI=sparse(u1n,u1n,chiI,n,n); 
    phiA=sparse(u1n,u1n,phiA_v,n,n); chiA=sparse(u1n,u1n,chiA,n,n); 
    
    % transmission terms 
    DeltaInv=sparse(u1n,u1n,1./(u*(St*MS+Et*ME+Pt*MP+It*MI+At*MA+Rt*MR)));
    %ThetaP=St*MS*(U-epsilon)*betaP*DeltaInv*MP';
    %ThetaI=St*MS*(U-epsilon)*betaI*DeltaInv*MI';
    %ThetaA=St*MS*(U-epsilon)*betaA*DeltaInv*MA';
    % NatComm beta assumption
    ThetaP=St*MS*(U-epsilon)*DeltaInv*MP'*betaP;
    ThetaI=St*MS*(U-epsilon)*DeltaInv*MI'*betaI;
    ThetaA=St*MS*(U-epsilon)*DeltaInv*MA'*betaA;
    
    % inverting transition rates for the next-generation matrix
    phiE_inv=sparse(u1n,u1n,1./phiE_v,n,n);
    phiP_inv=sparse(u1n,u1n,1./phiP_v,n,n);
    phiI_inv=sparse(u1n,u1n,1./phiI_v,n,n);
    phiA_inv=sparse(u1n,u1n,1./phiA_v,n,n);
    
    % next-generation-matrix
    K=deltaE*(ThetaP+sigma*deltaP*ThetaI*phiI_inv+(1-sigma)*deltaP*ThetaA*phiA_inv)*phiE_inv*phiP_inv;
    
    % effective reproduction number
    R_t=eigs(K,1,'largestabs');
    
    %  Jacobian matrix
    J0=[-mu*U Z -ThetaP -ThetaI -ThetaA Z Z Z Z Z;
        Z -phiE ThetaP ThetaI ThetaA Z Z Z Z Z; 
        Z deltaE*U -phiP Z Z Z Z Z Z Z;
        Z Z sigma*deltaP*U -phiI Z Z Z Z Z Z; 
        Z Z (1-sigma)*deltaP*U Z -phiA Z Z Z Z Z;
        Z chiE Z Z Z -phiEq*U Z Z Z Z;
        Z Z chiP Z Z deltaE*U -phiPq*U Z Z Z;
        Z Z Z eta*U+chiI Z Z sigma*deltaP*U -phiIh*U Z Z;
        Z Z Z Z chiA Z (1-sigma)*deltaP*U Z -phiAq*U Z; 
        Z Z Z gammaI*U gammaA*U Z Z gammaI*U gammaA*U -mu*U];
    % Hermitian matrix for epidemicity analysis
    H0_temp=W*J0*Wpi;
    H0=(H0_temp+H0_temp')/2;
    % effective epidemicity index
    e_t=eigs(H0,1,'largestreal');
end