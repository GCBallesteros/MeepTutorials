; Released under the Apache v2 License
; For a full copy of the license please visit:
; https://www.apache.org/licenses/LICENSE-2.0
;
; Author: Guillem Ballesteros
; Date: March 2018


; Gaussian Beam Function
(define (gaussian-beam zR k z) (lambda (x)
  (let
    ((z-prime (/ z zR))
     (k-norm (vector3-norm k))
     (r2 (vector3-dot x x)))
    (*
      (/ 1 (sqrt (+ 1 (* z-prime z-prime))))
      (exp (+
        (/ (* -1 k-norm r2)        (* 2 zR (+ 1 (* z-prime z-prime))))
        (/ (* 0-1i z-prime k-norm r2) (* 2 zR (+ 1 (* z-prime z-prime))))
        (* 0+1i (atan z-prime))
        (* 0-1i (vector3-dot k x))))))))