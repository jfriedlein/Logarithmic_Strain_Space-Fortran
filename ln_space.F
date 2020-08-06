c ###############################################################################
c ln-space - preprocessing 3D
c inputs:
c * hsv(NHV+1:NHV+9)
c outputs:
c * eigenvalues
c * eigenbasis
c * ea,da,fa
c * Hencky_strain_E 
      subroutine pre_ln( defoGrad_vec, eigenvalues, Ma,
     &                   ea, da, fa, Hencky_strain_E )
        use Tensor
        use TensorXLSDYNA
        implicit none
c Declarations        
        real, dimension(9), intent(in) :: defoGrad_vec
        type(Tensor2) :: Hencky_strain_E, F1
        type(Tensor2) :: right_CauchyGreen_C, dummy_C
        double precision eigenvectors_vec(3,3), eigenvalues(3), tmp_EW
        double precision ea(3), da(3), fa(3)
        type(Tensor2), dimension(3) :: Ma
        type(Tensor1), dimension(3) :: eigenvectors
        type(Tensor1) :: tmp_eigenvector
        integer a,b,c,i,j,k,l
c Extract the deformation gradient
        F1 = defoGrad( defoGrad_vec )
c @todo Check why printing the unsymmetric tensor F1 shows the transposed?
c Even though accessing the entry F1(1,3) contains the correct value
        right_CauchyGreen_C = transpose(F1) * F1
c Eigenvalues and eigenvectors of the right Cauchy-Green tensor
c The latter is always 3x3, hence we can use an optimized computation.
c @note Be aware that the function destroyes the inputted matrix hence
c we give a dummy variable, that could be deleted afterwards
        dummy_C = right_CauchyGreen_C
        call DSYEVJ3(dummy_C, eigenvectors_vec, eigenvalues)
c Store the eigenvectors as first order tensor to be able to handle them easier
c and e.g. apply the TTB dyadic product
        do a=1,3
          ! Be aware that the eigenvectors must be retrieved correctly
          eigenvectors(a) = eigenvectors_vec(1:3,a)
          enddo
c Sort the eigenvalues and eigenvectors in descendig order
c @todo-optimize: Find a better sort algo that also swapes the eigenvectors
c Negative EWs are not allowed, are they?. So we don't care about them
          if ( eigenvalues(1) > eigenvalues(2)
     &   .and. eigenvalues(2) > eigenvalues(3) ) then
              ! That's what we want
          ! 213 The second EW is larger than the first, but the third is still the smallest
          ! So we just swap the first and second
          elseif ( eigenvalues(2) > eigenvalues(1)
     &   .and. eigenvalues(1) > eigenvalues(3) ) then
              tmp_EW = eigenvalues(1)
              eigenvalues(1) = eigenvalues(2)
              eigenvalues(2) = tmp_EW
              tmp_eigenvector = eigenvectors(1)
              eigenvectors(1) = eigenvectors(2)
              eigenvectors(2) = tmp_eigenvector
          ! 321 The third EW is the largest and the second EW still the middle
          ! Swap the first and third one
          elseif ( eigenvalues(3) > eigenvalues(2)
     &   .and. eigenvalues(2) > eigenvalues(1) ) then
              tmp_EW = eigenvalues(1)
              eigenvalues(1) = eigenvalues(3)
              eigenvalues(3) = tmp_EW
              tmp_eigenvector = eigenvectors(1)
              eigenvectors(1) = eigenvectors(3)
              eigenvectors(3) = tmp_eigenvector  
          ! 132 The first EW is the first EW, but the second is the smallest
          ! Swap second and third EW
          elseif ( eigenvalues(1) > eigenvalues(3)
     &   .and. eigenvalues(3) > eigenvalues(2) ) then
              tmp_EW = eigenvalues(2)
              eigenvalues(2) = eigenvalues(3)
              eigenvalues(3) = tmp_EW
              tmp_eigenvector = eigenvectors(2)
              eigenvectors(2) = eigenvectors(3)
              eigenvectors(3) = tmp_eigenvector
          ! 231 The third is the largest, the second the smallest
          ! Store the third, move 1. and 2. back and insert the first
           elseif ( eigenvalues(3) > eigenvalues(1)
     &   .and. eigenvalues(1) > eigenvalues(2) ) then
              tmp_EW = eigenvalues(3)
              eigenvalues(3) = eigenvalues(2)
              eigenvalues(2) = eigenvalues(1)
              eigenvalues(1) = tmp_EW
              tmp_eigenvector = eigenvectors(3)
              eigenvectors(3) = eigenvectors(2)
              eigenvectors(2) = eigenvectors(1)
              eigenvectors(1) = tmp_eigenvector
          ! 312 The second is the largest, the first the smallest
          ! Store the second, move the third to the middle and the first to the end, insert the second
           elseif ( eigenvalues(2) > eigenvalues(3)
     &   .and. eigenvalues(3) > eigenvalues(1) ) then
              tmp_EW = eigenvalues(2)
              eigenvalues(2) = eigenvalues(3)
              eigenvalues(3) = eigenvalues(1)
              eigenvalues(1) = tmp_EW
              tmp_eigenvector = eigenvectors(2)
              eigenvectors(2) = eigenvectors(3)
              eigenvectors(3) = eigenvectors(1)
              eigenvectors(1) = tmp_eigenvector   
      endif

c @todo-extend Now we would have to do some things for 2D. So later on add these.
c @todo-optimize Find a better way to init the tensor (maybe automatially on creation)
        forall(i=1:3,j=1:3) Hencky_strain_E%ab(i,j)=0.
c
        do a=1,3
          ea(a) = 0.5 * log(eigenvalues(a))
          da(a) = 1. / eigenvalues(a)
          fa(a) = -2. * eigenvalues(a)**(-2)
c @todo Check whether the eigenvalues need to be sorted
          Ma(a) = (eigenvectors(a)) .dya. (eigenvectors(a))
          Hencky_strain_E = Hencky_strain_E + ea(a) * Ma(a)
      enddo
c        
!      if ( eps(2)>1e-12 ) then
!          write(*,*) "F1=", F1
!          write(*,*) "eigenvectors=", eigenvectors
!          write(*,*) "eigenvalues=", eigenvalues
!          write(*,*) "Hencky_strain=", Hencky_strain_E
!          pause
!      endif
c
        return
      end
      
      
c ###############################################################################
c ln-space - Compute projection tensor
c inputs:
c * eigenvalues
c * ea, da, fa
c * Ma
c output:
c * projection_tensor_P   
c * theta, xi, eta
      subroutine post_ln_projTensor( eigenvalues, Ma,
     &            ea, da, fa, projection_tensor_P, theta, xi, eta )
        use Tensor
        implicit none
c Declarations        
        double precision eigenvalues(3)
        double precision ea(3), da(3), fa(3)
        type(Tensor2), dimension(3), intent(in) :: Ma
        type(Tensor4) :: projection_tensor_P, tensor_operator_G, tmp
        type(Tensor2) :: theta, xi
        double precision eta
        double precision comp_tolerance
        integer a,b,c,i,j,k,l        
c USER-Parameters
        comp_tolerance = 1e-8
c Set up the coefficients theta, xi and eta
       ! For three different eigenvalues \f$ \lambda_a \neq \lambda_b \neq \lambda_c \f$
	 if (
     &     (.NOT.(abs(eigenvalues(1)-eigenvalues(2)) < comp_tolerance) )
     &     .AND.
     &     (.NOT.(abs(eigenvalues(1)-eigenvalues(3)) < comp_tolerance) )
     &     .AND.
     &     (.NOT.(abs(eigenvalues(2)-eigenvalues(3)) < comp_tolerance) )
     &	 ) then
	    eta = 0.0
          do a=1,3
              do b=1,3
 				if (a /= b) then
					theta%ab(a,b) = (ea(a) - ea(b))
     &  							  / (eigenvalues(a) - eigenvalues(b))
					xi%ab(a,b) = (theta%ab(a,b) - 0.5 * da(b))
     &							   / (eigenvalues(a) - eigenvalues(b))

                     do c=1,3
						if ((c /= a) .AND. (c /= b)) then
							eta = eta +
     &									ea(a)
     &									/ (2.0
     &									   * (eigenvalues(a)
     &										  - eigenvalues(b))
     &									   * (eigenvalues(a)
     &										  - eigenvalues(c)))
                          endif
                     enddo
                 endif
              enddo
          enddo
      ! For three equal eigenvalues \f$ \lambda_a = \lambda_b = \lambda_c \f$    
      elseif ( (abs(eigenvalues(1) - eigenvalues(2)) < comp_tolerance)
     &          .AND.
     &          (abs(eigenvalues(2)-eigenvalues(3)) < comp_tolerance) )
     & then
		eta = 0.0
          do a=1,3
              do b=1,3
				if (a /= b) then
					theta%ab(a,b) = 0.5 * da(1)
					xi%ab(a,b) = (1.0 / 8.0) * fa(1)
				endif
              enddo
          enddo
		eta = (1.0 / 8.0) * fa(1)
      ! For two equal eigenvalues a and b: \f$ \lambda_a = \lambda_b \neq \lambda_c \f$
      elseif ( (abs(eigenvalues(1) - eigenvalues(2)) < comp_tolerance)
     &		 .AND.
     &		 ( .NOT.(abs(eigenvalues(2)-eigenvalues(3)) 
     &                < comp_tolerance) ) )
     & then
		eta = 0.0
          do a=1,3
              do b=1,3
				if ((a /= b) .AND. ((a == 3) .OR. (b == 3))) then
					theta%ab(a,b) = (ea(a) - ea(b))
     &							  / (eigenvalues(a) - eigenvalues(b))
					xi%ab(a,b) = (theta%ab(a,b) - 0.5 * da(b))
     &						     / (eigenvalues(a) - eigenvalues(b))
				endif
              enddo
          enddo
		theta%ab(1,2) = 0.5 * da(1)
		theta%ab(2,1) = theta%ab(1,2)
		xi%ab(1,2) = (1.0 / 8.0) * fa(1)
		xi%ab(2,1) = xi%ab(1,2)
		eta = xi%ab(3,1)
      ! For two equal eigenvalues a and c: \f$ \lambda_a = \lambda_c \neq \lambda_b \f$
      elseif ( (abs(eigenvalues(1) - eigenvalues(3)) < comp_tolerance)
     &		  .AND.
     &		  (.NOT.(abs(eigenvalues(2)-eigenvalues(3)) < comp_tolerance)) )
     & then
          eta = 0.0
          do a=1,3
              do b=1,3
				if ( (a /= b) .AND. ((a == 2) .OR. (b == 2)) ) then
					theta%ab(a,b) = (ea(a) - ea(b))
     &							  / (eigenvalues(a) - eigenvalues(b))
					xi%ab(a,b) = (theta%ab(a,b) - 0.5 * da(b))
     &						   / (eigenvalues(a) - eigenvalues(b))
				endif
              enddo
          enddo
		theta%ab(1,3) = 0.5 * da(1)
		theta%ab(3,1) = theta%ab(1,3)
		xi%ab(1,3) = (1.0 / 8.0) * fa(1)
		xi%ab(3,1) = xi%ab(1,3)
		eta = xi%ab(2,1)
	! For two equal eigenvalues b and c: \f$ \lambda_b = \lambda_c \neq \lambda_a \f$
      else if ( (abs(eigenvalues(2) - eigenvalues(3)) < comp_tolerance)
     &		  .AND.
     &		  (.NOT.(abs(eigenvalues(1)-eigenvalues(2)) < comp_tolerance)))
     & then
		eta = 0.0
          do a=1,3
              do b=1,3
				if ( (a /= b) .AND. ((a == 1) .OR. (b == 1)) ) then
					theta%ab(a,b) = (ea(a) - ea(b))
     &							  / (eigenvalues(a) - eigenvalues(b))
					xi%ab(a,b) = (theta%ab(a,b) - 0.5 * da(b))
     &						   / (eigenvalues(a) - eigenvalues(b))
				endif
              enddo
          enddo
		theta%ab(2,3) = 0.5 * da(2)
		theta%ab(3,2) = theta%ab(2,3)
		xi%ab(2,3) = (1.0 / 8.0) * fa(2)
		xi%ab(3,2) = xi%ab(2,3)
		eta = xi%ab(1,2)
	 else
          write(*,*) "Found a non distinguisable eigenvalue triple"
          write(*,*) "eigenvalues=", eigenvalues
          pause
       endif
c
c Set the fourth order projection tensor to transform the stress tensor
      forall(i=1:3, j=1:3, k=1:3, l=1:3)
     &      projection_tensor_P%abcd(i,j,k,l) = 0.
       do a=1,3
c The variable 'tmp' seems to be necessary, else we cannot assign the value of tmp
          tmp = da(a) * Ma(a).dya.Ma(a)
          projection_tensor_P = projection_tensor_P + tmp
          do b=1,3
			if (b /= a) then
                  forall (i=1:3, j=1:3, k=1:3, l=1:3)
     &                tensor_operator_G%abcd(i,j,k,l) =
     &                                    Ma(a)%ab(i,k) * Ma(b)%ab(j,l)
     &                                    +Ma(a)%ab(i,l) * Ma(b)%ab(j,k)
				projection_tensor_P = projection_tensor_P +
     &                        theta%ab(a,b) * tensor_operator_G
              endif
	    enddo
      enddo
c        
      if ( .false.) then
          write(*,*) "projection_tensor_P=",projection_tensor_P
          write(*,*) "theta=",theta
          write(*,*) "xi=",xi
          write(*,*) "eta=",eta
      endif
        return
      end

      
   
c ###############################################################################
c ln-space - utan postprocessing
c inputs:
c * tangent_C
c * stress_S
c * eigenvalues
c * ea, da, fa
c * Ma
c output:
c * tangent_E
c @todo ensure that we get stress_S in utan when calling sig()
      subroutine post_ln_utan( tangent_C, stress_T, Ma, 
     &   projection_tensor_P, theta, xi, eta, fa, tangent_E )
        use Tensor
        implicit none
c Declarations        
        type(Tensor4), intent(in) :: tangent_C
        type(Tensor2), intent(in) :: stress_T
        type(Tensor2), dimension(3), intent(in) :: Ma
        double precision, dimension(3), intent(in) :: fa
        type(Tensor4) :: projection_tensor_P, T_DC_L, tmp
        type(Tensor4) :: tangent_E
        type(Tensor4) :: ten_r_abb, ten_l_abb, ten_r_bba,
     &              ten_l_bba, ten_r_bab, ten_l_bab,
     &                ten_r_abc, ten_l_bca
        integer :: i, j, k, l, a, b, c
        type(Tensor2) :: theta, xi, Ma_a, Ma_b, Ma_c
        type(Tensor2) :: aaa, bbb, ccc
        double precision :: eta
c        
c Compute T:L
        forall(i=1:3, j=1:3, k=1:3, l=1:3)
     &      T_DC_L%abcd(i,j,k,l) = 0.
c      
       do a=1,3
		tmp = fa(a)
     &				* (stress_T ** Ma(a))
     &				* (Ma(a).dya.Ma(a))
          T_DC_L = T_DC_L + tmp
c
          do b=1,3
			if (b /= a) then
                  call ten_op_F_right( Ma, a, b, b,
     &             stress_T, ten_r_abb )
                  call ten_op_F_left( Ma, a, b, b,
     &             stress_T, ten_l_abb )
                  call ten_op_F_right( Ma, b, b, a,
     &             stress_T, ten_r_bba )
                  call ten_op_F_left( Ma, b, b, a,
     &             stress_T, ten_l_bba )
                  call ten_op_F_right( Ma, b, a, b,
     &             stress_T, ten_r_bab )
                  call ten_op_F_left( Ma, b, a, b,
     &             stress_T, ten_l_bab )
c
				tmp = 2.0 * xi%ab(a,b)
     &					 * ( ten_r_abb + ten_l_abb + ten_r_bba
     &                         + ten_l_bba + ten_r_bab + ten_l_bab )
c                  
                  T_DC_L = T_DC_L + tmp
				do c=1,3
                      Ma_c = Ma(c)
					if ( (c /= a) .AND. (c /= b) ) then
                          call ten_op_F_right( Ma, a, b, c,
     &                        stress_T, ten_r_abc )
                          call ten_op_F_left( Ma, b, c, a,
     &                        stress_T, ten_l_bca )
						tmp = 2.0 * eta * ( ten_r_abc + ten_l_bca )
                          T_DC_L = T_DC_L + tmp
					endif
                  enddo                     
			endif
          enddo
       enddo
c        
        tangent_E = projection_tensor_P**tangent_C**projection_tensor_P
     &                + T_DC_L
c       
        return
              end
              
c ###############################################################################
      subroutine ten_op_F_right( Mabc, a, b, c, stress_T, ten_x_xxx )
          use Tensor
          implicit none

          type(Tensor2), dimension(3) :: Mabc
          type(Tensor4) :: ten_x_xxx
          type(Tensor2) :: temp_tensor, MbTMc, stress_T
          integer i, j, k, l, a, b, c
c          
	   temp_tensor = stress_t * Mabc(c)
	   MbTMc = Mabc(b) * temp_tensor

        forall(i=1:3,j=1:3,k=1:3,l=1:3)
     &	 ten_x_xxx%abcd(i,j,k,l) = 
     &       Mabc(a)%ab(i,k)*MbTMc%ab(j,l)+Mabc(a)%ab(i,l)*MbTMc%ab(j,k)
        
      return
      end
        
c ###############################################################################
      subroutine ten_op_F_left( Mabc, a, b, c, stress_T, ten_x_xxx )
          use Tensor
          implicit none

          type(Tensor2), dimension(3) :: Mabc
          type(Tensor4) :: ten_x_xxx
          type(Tensor2) :: temp_tensor, MaTMb, stress_t
          integer i, j, k, l, a, b, c
c          
	   temp_tensor = stress_t * Mabc(b)
	   MaTMb = Mabc(a) * temp_tensor

        forall(i=1:3,j=1:3,k=1:3,l=1:3)
     &	 ten_x_xxx%abcd(i,j,k,l) = 
     &       MaTMb%ab(j,l)*Mabc(c)%ab(i,k)+MaTMb%ab(j,k)*Mabc(c)%ab(i,l)
c        
      return
      end
              
c ###############################################################################
c ln-space - postprocessing
c inputs:
c * stress_T
c * tangent_C
c * eigenvalues
c * ea, da, fa
c * Ma
c output:
c * stress_S
c * tangent_E
       subroutine post_ln( stress_T, tangent_C, eigenvalues, Ma,
     &              ea, da, fa, stress_S, tangent_E )     
        use Tensor
        implicit none
c Declarations        
        type(Tensor2), intent(in) :: stress_T
        double precision, intent(in) :: eigenvalues(3)
        double precision, intent(in) ::  ea(3), da(3), fa(3)
        type(Tensor2), dimension(3), intent(in) ::  Ma
        type(Tensor2) :: stress_S
        type(Tensor4) :: projection_tensor_P, tangent_C, tangent_E
        type(Tensor2) :: theta, xi
        double precision eta
c
        call post_ln_projTensor ( eigenvalues, Ma,
     &            ea, da, fa, projection_tensor_P, theta, xi, eta )
c Transform the stress
c The double contraction is performed by the '**' operator
      stress_S = stress_T ** projection_tensor_P
      
c Transform the tangent
      call post_ln_utan( tangent_C, stress_T, Ma, 
     &   projection_tensor_P, theta, xi, eta, fa, tangent_E )     
c
      return
      end
              