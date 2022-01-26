 module Propagation_m
 USE NumParameters_m
 USE Op_m


  SUBROUTINE TAYLOR()
   USE NumParameters_m
    USE Basis_m 
   USE psi_m, ONLY : psi_t 
  
   
    !TYPE(Besis_m):: Besis
    TYPE (psi_t)           :: psi, ans
    complex (kind=Rk)   :: vec(6),vec0(6)
    real(kind= Rk)      :: t(12) 
    integer:: k, n, kk, Ntps  !(Ntps = nombre de points de grille temporelle)
    real (kind=Rk),parameter:: sigma, scaletps,PI,k0
          real(kind= Rk) :: Rkk,NN
      !allocate(vec(size(Besis%x),vec0(size(Besis%x))
             !allocate(t(N))
             k0 = ONE ;  sigma = FOURTH ; scaletps= ONETENTH; Ntps = 12
             
             
             ! scaletps = pas temporel
             ! sigma = largeur de la gaussienne
             !k0 = p0 vitesse initiale
                        !++++++++++++++++++++initialisation+++++++++++++
          
vec0(:) = (ONE/(PI*sigma*sigma)exp(-EYE*k0*Besis%x)
  vec0(:) vec0(:)*exp(-((Besis%x-ZERO)*(Besis%x-ZERO))/2*sigma*sigma)
     t(:) = (n-1)*scaletps
     ans(:) = vec0(:)
     psi(:) = ZERO
      kk = ONE
       Rkk = ONE
     
     Do n = 1,Ntps,1
        Do k = 1, 1000
     Rkk = ONE
      vec(:) = vec0(:)
          CALL  calc_OpPsi(Op,vec,vec)
       vec(:) = -EYE*vec(:)
        kk = kk+1 .AND. Rkk = Rkk*(scaletps/kk)
        ans(:) = ans(:) + Rkk*vec(:)
        NN= Rkk*scaleQ*sqrt(matmul(vec(:),transpose(vec(:))))
         IF(NN .lt. TEN**(-TEN ) ) then
             exit
          enddo   
     ENDO     
        
        
       
      
    
   END SUBROUTINE TAYLOR  
   
   
   
   end module Propagation_m
