; Modified Meep ring resonator example to show how to make
; a frequency monitor. I added an output port so that energy
; dissipated out of the system faster.
;
; Original code borrowed from:
; https://github.com/stevengj/meep/blob/master/examples/ring.ctl

; Calculating 2d ring-resonator modes, from the Meep tutorial.

(define-param n 3.4) ; index of waveguide
(define-param w 1) ; width of waveguide
(define-param r 1) ; inner radius of ring
(define-param output-gap 0.1)

(define-param pad 5) ; padding between waveguide and edge of PML
(define-param dpml 2) ; thickness of PML

(define nfreq 200)

(define sxy (* 2 (+ r w pad dpml))) ; cell size
(set! geometry-lattice (make lattice (size sxy sxy no-size)))

; Create a ring waveguide by two overlapping cylinders - later objects
; take precedence over earlier objects, so we put the outer cylinder first.
; and the inner (air) cylinder second. Add an output port.
(set! geometry (list
		(make cylinder (center 0 0) (height infinity)
		      (radius (+ r w)) (material (make dielectric (index n))))
		(make cylinder (center 0 0) (height infinity)
		      (radius r) (material air))
		(make block (center 0 (+ r w (/ w 2) output-gap))
			  (size infinity w) (material (make dielectric (index n))))))

(set! pml-layers (list (make pml (thickness dpml))))
(set-param! resolution 10)

; If we don't want to excite a specific mode symmetry, we can just
; put a single point source at some arbitrary place, pointing in some
; arbitrary direction.  We will only look for TM modes (E out of the plane).

(define-param fcen 0.25) ; pulse center frequency
(define-param df 0.15)  ; pulse width (in frequency)
(set! sources (list
               (make source
                 (src (make gaussian-src (frequency fcen) (fwidth df)))
                 (component Ez) (center (+ r 0.1) 0))))

;; = = = = = = = = = = = = = = = = = = = = = = = = = =
; Add metadata to H5 file helpers
; HDF5 METADATA Utils
(define (attr->string attr-pair)
  (string-append " -attr " (car attr-pair) " " (number->string (exact->inexact (second attr-pair)))))

(define (add-attributes fname attributes)
  (let
    ((command (string-append "python ./add_attr.py -fname " (get-filename-prefix) "-" fname ".h5 "))
      (attribute-string (string-join (map attr->string attributes))))
    (string-append command attribute-string)))

; Helpers to create monitors
(define (make-volume x-min x-max y-min y-max)
  (volume
    (center (/ (+ x-max x-min) 2) (/ (+ y-max y-min) 2) 0)
    (size (- x-max x-min) (- y-max y-min) 0)))

(define (freq-field-monitor monitor-name dt vol output)
  (at-every dt (to-appended monitor-name (in-volume vol output))))
;; = = = = = = = = = = = = = = = = = = = = = = = = = =

; Create the monitor and the step
(define dt-monitor (let ((fmax (+ fcen (/ df 2))))
        (/ 1 (* 2 fmax))))
(define ring-volume (make-volume (* -1 (+ r w 0.2)) (+ r w 0.2)
                                 (* -1 (+ r w 0.2)) (+ r w 0.2)))
(define ring-monitor (freq-field-monitor "freq-monitor" dt-monitor ring-volume output-efield-z))


; Run simulation
(run-sources+ (stop-when-fields-decayed 25 Ez (vector3 (+ r 0.1) 0) 1e-8)
	      (at-beginning output-epsilon)
	      ring-monitor
	      (dft-ldos fcen df nfreq))

; Add metadata to output
(system (add-attributes "freq-monitor" (list (list "xmin" (* -1 (+ r w 0.2)))
                                             (list "xmax" (+ r w 0.2))
                                             (list "ymin" (* -1 (+ r w 0.2)))
                                             (list "ymax" (+ r w 0.2))
                                             (list "dt" dt-monitor)
                                             (list "minfreq" (- fcen (/ df 2)))
                                             (list "maxfreq" (+ fcen (/ df 2)))
                                             (list "res" 10))))

; Meep gets stuck in interactive mode after a system call. Don't know how to fix
; it yet so I just force quit.
(quit)
