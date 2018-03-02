; Released under the Apache v2 License
; For a full copy of the license please visit:
; https://www.apache.org/licenses/LICENSE-2.0
;
; Author: Guillem Ballesteros
; Date: February 2018

; Include definition of the source
(include "./gaussian_source.scm")

; Bounding Box
(define-param s 20)
(define-param dpml 1)
(define-param T 50) ; run time
(define sxy (+ s (* 2 dpml))) ; cell size, including PML

(set! geometry-lattice (make lattice (size sxy sxy no-size)))
(set! pml-layers (list (make pml (thickness dpml))))
(set-param! resolution 10)

; Source Parameters
(define-param fcen 0.8)
(define-param df 0.001)

; Beam Parameters
(define-param theta (deg->rad 0)) ; Beam direction positive is downwards
(define-param w0 1.5) ; Beam waist
(define-param n 1) ; Index of refraction of medium
(define-param dist-to-waist 0) ; positive focuses on front of the source

(define zR (* pi n w0 w0 fcen)) ; Rayleigh range derived from beam waist

(define k (vector3-scale
  (* 2 pi fcen)
  (unit-vector3 (vector3 (cos theta) (sin theta)))))

(set! sources
      (list
        (make source
      	  (src (make continuous-src (frequency fcen) (fwidth df)))
      	  (component Ez)
          (center (* -0.48 s) 0)
          (size 0 s)
      	  (amp-func (gaussian-beam zR k dist-to-waist)))))

(run-until T (at-end (output-png Ez "-S3 -Zc bluered")))
