(rt:deftest s*m
    (let ((B #A((2 3) single-float
		(1.0 0.4 0.2)
		(0.1 0.2 0.3))))
      (values
       (with-matrixes (* 12.0 B) )
       (with-matrixes (* 12.0 B))))
  #1=#A((2 3) single-float (12.0 4.8 2.4) (1.2 2.4 3.6000001))
  #1#)
(regression-test:deftest scalar
    (with-matrixes 12) 12)

(rt:deftest scalar*
    (with-matrixes (* 12)) 12)

(rt:deftest scalar*3
    (with-matrixes (* 12 12)) 144)


(rt:deftest trace*
    (let ((B (make-array '(5 7) :element-type 'single-float :initial-element 3.0))
	  (C (make-array '(7 5) :element-type 'single-float :initial-element 1.0))
	  (alpha 0.5))
      (with-matrixes
	  (trace (* alpha B C))))
  52.5)

(rt:deftest mult
    (let ((B (make-array '(5 7) :element-type 'single-float :initial-element 3.0)))
      (with-matrixes
	  (* B (transpose B))))
  #.(make-array '(5 5) :element-type 'single-float :initial-element 63.0))

(rt:deftest mult2
    (let ((B #A ((3 2) single-float
		 (1.0 0.0)
		 (0.0 0.2)
		 (0.1 2.5))))
      (with-matrixes
	  (* B (diag 1 2))))
  #2A((1.0 0.0) (0.0 0.4) (0.1 5.0)))

;; rewrite of tests in multiply-test
(rt:deftest direct2
    (with-matrixes
	(*
	 #2A((2 1 4 5)(0 -3 -1 7) (6 2 9 -8))
	 #2a((3  6)(-1 1) (5 0) (2 -4)))
	:field fixnum)
  #2A((35 -7) (12 -31) (45 70)))

(rt:deftest transposed2
    (with-matrixes
	(* (transpose #2A((2 0 6)(1 -3 2)(4 -1 9)(5 7 -8)))
	   #2a((3 6)(-1 1)(5 0)(2 -4)) )
	:field fixnum)

  #2A((35 -7)(12 -31)(45 70)))

(with-matrixes
    (+ 1.0 2.0 3.0))

(rt:deftest (linear-combination mat)
    (with-matrixes
	(+ (* 3 #.(make-array '(3 3) :initial-element 1))
	   (* 2 #2a((1 2 3)
		    (4 5 6)
		    (7 8 9))))
	:field fixnum)
  #2A ((5 7 9) (11 13 15) (17 19 21)))

(rt:deftest (linear-combination-mismatch mat)
  (handler-case
      (with-matrixes
	  (+ (* 1 #2A ((1 2 3)))
	     (* 1 #2a ((1 2))))
	  :field fixnum)
    (t (e) (type-of e)))
  simple-error)