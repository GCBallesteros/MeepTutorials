(define-param structure? true)

; Periodic Boundary Conditions
(set-param! k-point (vector3 0 0 0))

; Source
(define-param min-freq 0.5);0.95
(define-param max-freq 1.5);1.12
(define fcen (/ (+ min-freq max-freq) 2))
(define df (- max-freq min-freq))

; DBR
(define-param n1 3)
(define-param n2 2)
(define-param n-pairs 20)
(define-param central-wl 1)

; Total thickness of DBR
(define (thickness-DBR n1 n2 n-pairs wl)
  (*
    n-pairs
    wl
    0.25
    (+ (/ 1 n1) (/ 1 n2))))

; Return list of blocks that makes the DBR
(define (add-DBR n1 n2 n-pairs wl y-pos)
  (let
    ((layer1-t (* 0.25 (/ wl n1)))
     (layer2-t (* 0.25 (/ wl n2)))
     (make-layer (lambda (layer-t index top) (make block
                                                  (center 0 (- top (/ layer-t 2)) 0)
                                                  (material (make dielectric (epsilon (* index index))))
                                                  (size infinity layer-t infinity)))))
    (let
      ((make-pair (lambda (top) (list
              (make-layer layer1-t n1 top)
              (make-layer layer2-t n2 (- top layer1-t))))))
      (do
        ((ii 0 (+ ii 1))
         (top y-pos (- top layer1-t layer2-t))
         (layer-list  '() (append layer-list (make-pair top))))
        ((= ii n-pairs) layer-list)))))

; Lattice
; High res because convergence to analytical answer is slowww
(set-param! resolution 60)
(define sx 4)
(define sy (+
  (thickness-DBR n1 n2 n-pairs central-wl)
  6))
(define dpml 1)

(set! geometry-lattice (make lattice (size sx (+ sy (* 2 dpml)) no-size)))
(set! pml-layers (list (make pml (thickness dpml) (direction Y))))

(define top-mirror (/ (thickness-DBR n1 n2 n-pairs central-wl) 2))

(set! geometry 
    (if structure?
      (append 
        (add-DBR n1 n2 n-pairs central-wl top-mirror)
        '())
      '()))

; Sources
;Plane Wave Source functoin
(define (plane-wave k) (lambda (x)
  (exp (* 0-1i (vector3-dot k x)))))

;Propagation Constant
(define-param theta (deg->rad 90)) ; Beam direction positive is downwards
(define k (vector3-scale
  (* 2 pi fcen)
  (unit-vector3 (vector3 (cos theta) (sin theta)))))

(set! sources
      (list
        (make source
          (src (make gaussian-src (frequency fcen) (fwidth df)))
          (component Ez)
          (center 0 (/ sy 2))
          (size sx 0)
          (amp-func (plane-wave k)))))

; Set transmission monitor
(define nfreq 200)
(define top                                               
      (add-flux fcen df nfreq
              (make flux-region
              (center 0 (+ top-mirror 0.5) 0)
              (size sx 0 sx))))

(define bottom                                               
      (add-flux fcen df nfreq
              (make flux-region
              (center 0 (* -1 (+ top-mirror 0.5)) 0)
              (size sx 0 sx))))

(if structure?
  (load-minus-flux "ref_flux" top))

(run-sources+ (stop-when-fields-decayed 25 Ez (vector3 0 (/ sy 2) 0) 1e-8)
              (at-beginning output-epsilon))

(if (not structure?)
  (save-flux "ref_flux" top))

(display-fluxes top bottom)
