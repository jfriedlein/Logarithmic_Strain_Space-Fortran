# Logarithmic_Strain_Space-Fortran
Fortran code for the transformation into the logarithmic strain space (ln-space)

(also available for C++ library deal.II [here](https://github.com/jfriedlein/Logarithmic_Strain_Space-dealii))

## Requirements
* The here published version of the ln-space is implemented in tensor notation using the tensor toolbox ttb that is also required and available [here](https://github.com/adtzlr/ttb).
* Moreover, for convenience the ttb extension for LS-Dyna (ttbXLSDYNA) needs to be used, which, at the moment, is only available [here](https://github.com/jfriedlein/ttb/tree/extension_test_XLSDYNA/ttbXLSDYNA). However, this is only being used to transform the deformation gradient from the LS-Dyna list storage into a second order tensor, which is trivial if you know the storage convention (see LS-Dyna user manual, appendix A).
* The example below uses some more modules, namely the [history hsv-manager](https://github.com/jfriedlein/history_hsv-manager_LS-Dyna), [material parameter cm-manager](https://github.com/jfriedlein/material-parameter_manager_LS-Dyna) and the ttbXkinematics extension for the stress and tangent push forwards. The former modules simplify the access of history variables and material parameters.


## Note
This code is not optimised. We are probably able to speed things up in many places (recompute/save data, Voigt notation, ...). The ln-space can also be set up such that it directly outputs the spatial stress and tangent. Unfortunately, the current implementation does not enable a simple switch between material and spatial output, because it is based on an older ln-space paper. So, it is simpler to do a separate push-forward in the end.

## The goal/When to use this code
The logarithmic strain space (herein often abbreviated as ln-space) is a very simple (in terms of the application) way to apply small strain material models to finite strains. So, in case you have a small strain model that you would like to use for applications exposed to large deformations/finite strains, the ln-space it probably the easiest way to achieve this.

All you need are three steps as schematically (simplified) shown in the figure:
1. Pre-processing from the world of finite strains into the logarithmic strain space
2. Calling your small strain model with the Hencky strain computed in the pre-processing
3. Post-processing the computed stress and tangent(s) by transforming them from the ln-space back into the real world.
And the best is, steps 1 and 3 are always the same (copy-paste) and the second step is the same as done for the small strain implementation.

<img src="https://github.com/jfriedlein/Logarithmic_Strain_Space-Fortran/blob/master/images/ln-space%20-%20overview.png" width="700">

Drawbacks?

Especially step 3, the post-processing is expensive independent of your material model. So, as far as efficiency is concerned, a simple material model utilising the ln-space is most certainly slower, than its derived finite strain equivalent (a model that was developed for finite strains). However, for complicated material models it can be faster to use the small strain model in the ln-space, instead of a finite strain equivalent (and it requires no further derivations/development to get from small to finite strains). Another disadvantage, is the limitation to small elastic strains. The latter is usually satisfied for metal forming and similar applications, where elastic strain are small and large plastic strain occur.

## Background
@todo add some equations

The transformation consists of three steps. First, we transform the deformation gradient `F` into the logarithmic strain space (ln-space) and obtain the Hencky strain (preprocessing). Secondly, the standard small strain material model is called using the Hencky strain and the usual history. The outcome is the stress measure `T` and the fourth order tangent `C`. Finally, we transform the stress and tangent modulus back from the logarithmic strain space to obtain e.g. the Second Piola-Kirchhoff stress tensor `S` and the Lagrangian tangent modulus `L` (postprocessing).

## Interface/How to use this code
To show how easy the usage is and to give you the boilerplate code snippets, here an example from an LS-DYNA UMAT (For details on LS-Dyna umats see [here](https://github.com/jfriedlein/usrmat_LS-Dyna_Fortran)).

We implement our user-defined material model in the file "dyn21umats.F", were we start with the includes. This ln-space is based on tensors, so we use the [ttb tensor toolbox for Fortran](https://github.com/adtzlr/ttb).

[dyn21umats.F]
```fortran
c Tensor Toolbox for Fortran by Andreas Dutzler
#define NOR4
#include '../ttb/ttb/ttb_library.F'
```
As mentioned above, we require the eigenproblem solver developed by Joachim Kopp, optimised for symmetric 3x3 matrices.
```fortran
c Subroutine for the computation of the eigenproblem by Joachim Kopp
c [https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html]
      include '../dsyevj3-F-1.0/dsyevj3.f'
```
And of course we have to include the ln-space itself, as follows.
```fortran
c Logarthmic strain space
      include "../Logarithmic_Strain_Space-Fortran/ln_space.F"
```

In the desired umat (here umat44), we follow the three steps (preprocessing, small strain model, postprocessing).
Everything until the next text block is boilerplate code, so can be reused without changes.

```fortran
      subroutine umat44 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,nnpcrv,cma,qmat,elsiz,idele,reject)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c Elastoplasticity - finite strain
c******************************************************************
c Thanks to the tensor module we can utilise our tensorial equations, but LS-Dyna
c still only provides (input) and desires (output) vectorial quantities.
      use Tensor
      use TensorXkinematics
      use hsv_manager
      use cm_manager
c
      include 'nlqparm'
      include 'bk06.inc'
      include 'iounits.inc'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
      integer nnpcrv(*)
      character*5 etype
      logical failel,reject
      INTEGER8 idele
c
c Declarations:
c Arrays must be declared as such
c But nevertheless declare every single variable,
c else you will get utter bs at some point.
      type(Tensor2) :: Hencky_strain_E,
     &                 stress_T, stress_S, stress_Cauchy,
     &                 defoGrad_F
      type(Tensor4) :: tangent_C, tangent_E, tangent_EE
      double precision eigenvalues(3)
      double precision ea(3), da(3), fa(3)
      type(Tensor2), dimension(3) :: Ma
      type(Tensor2) :: Eye
c
      double precision, dimension(6,6) :: es
      real :: fstrain
      logical failedIt
      failedIt = .false.
c      
c ###############################################################################
c ln-space - preprocessing
c inputs:
c * hsv ! contains the deformation gradient via flag IHYPER=1
c outputs:
c * eigenvalues
c * eigenbasis 'Ma'
c * ea, da, fa
c * Hencky_strain_E
c
c Finite strain (via ln-space) on/off
      fstrain = cm_get('finiteStrain____',cm)
c
      if ( fstrain>0.5 ) then
      call pre_ln( hsv_get_unsymTen2('defoG',hsv), 
     &             eigenvalues, Ma,ea, da, fa, Hencky_strain_E )
      else ! small strain
          Eye = identity2(Eye)
          Hencky_strain_E =
     &             symmetrize(hsv_get_unsymTen2('defoG',hsv) - Eye)
      endif
c
c ###############################################################################
```

Now step 2 begins, where we call our small strain material model with the computed Hencky_strain_E. In this example, the actual material model was outsourced into a separate file and hence the subroutine is called via "call".

```fortran
c Be aware that we now call our small strain material model with
c the "pseudo" small strain Hencky strain.
c The latter is, however, 
c a total strain, because it is computed from the deformation gradient
c instead of the incremental strain 'eps'. As a consequence the equations
c will differ from an LS-Dyna small strain model
c
c ###############################################################################
c Small strain material model
c inputs:
c * material parameters
c * Hencky_strain_E
c * history
c outputs:
c * stress
c * tangent
c * history
c
      call elasto_plasticity_code ( Hencky_strain_E, hsv, cm,
     &     stress_T, tangent_C, failedIt )
      
c
c ###############################################################################
c
c This ends our actual small strain material model, but before we can take
c the resulting stress tensor and go home, we have to transform it from the 
c logarithmic strain space (ln-space) back into the real world (post). Which is
c admittedly a complicated and expensive task.
```

And boilerplate code again, where we postprocess the stress and tangent and push
them forward, because LS-Dyna wants the spatial quantities.

```fortran
c ###############################################################################
c ln-space - postprocessing
c inputs:
c * eigenvalues
c * ea, da, fa
c * Ma
c output:
c * stress_S 
c * tangent_E 
c
      if ( fstrain>0.5 ) then
      call post_ln( stress_T, tangent_C, eigenvalues, Ma,
     &              ea, da, fa, stress_S, tangent_E )
      else
          stress_S = stress_T
          tangent_E = tangent_C
      endif
c
c ###############################################################################
c
c Save the stress into the return argument      
      sig(1:6) = asarray(voigt(stress_S),6)
c
c Push-Forward
       if ( fstrain>0.5 ) then          
         defoGrad_F = hsv_get_unsymTen2('defoG',hsv)

         ! Push-Forward of the 2. Piola-Kirchhoff stress
         stress_Cauchy = stress_push_forward( stress_S,
     &                                      defoGrad_F )
         sig(1:6) =  asarray(voigt(stress_Cauchy),6)
       endif
c      
c Push-Forward of the tangent
      es(1:6,1:6) = asarray(voigt(tangent_E),6,6)
c         
      if ( fstrain>0.5 ) then
          tangent_EE = tangent_push_forward( tangent_E,
     &                                       defoGrad_F )
          es(1:6,1:6) = asarray(voigt(tangent_EE),6,6)
      endif
c            
      call hsv_set_6x6(es,'es___',hsv)
c
      return
      end
```


## ToDo:
* test the use of include dsyevj in ln_space.F
* remove the need for all these modules, to make the code and the example standalone

Add a note (where to save it) and the download link
* "always existing lines of code" (maybe colored blocks)
