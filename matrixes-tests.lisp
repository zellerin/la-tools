(in-package linear-algebra)
(rt:deftest s*m
    (let ((B (with-matrixes (:constant (1 0.4 0.2 0.1 0.2 0.25) 2 3))))
      (values
       (with-matrixes (* 12.0 B) )))
  #.(with-matrixes (:constant (12 4.8 2.4 1.2 2.4 3) 2 3)))

(rt:deftest scalar
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

#+nil(rt:deftest mult2
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
	 #2a((3 6)(-1 1) (5 0) (2 -4)))
	:field fixnum)
  #2A((35 -7) (12 -31) (45 70)))

(rt:deftest transposed2
    (with-matrixes
	(* (transpose
	    #2A((2 0 6)(1 -3 2)(4 -1 9)(5 7 -8)))
	   #2a((3 6)(-1 1)(5 0)(2 -4)) )
	:field fixnum)

  #2A((35 -7)(12 -31)(45 70)))

(rt:deftest (linear-combination mat)
    (with-matrixes
	(+ (* 3 #.(make-array '(3 3) :initial-element 1))
	   (* 2 #2a((1 2 3)
		    (4 5 6)
		    (7 8 9))))
	:field fixnum)
  #2A ((5 7 9) (11 13 15) (17 19 21)))

(rt:deftest map
    (let ((a (make-array '(3 3) :initial-element 1))
	  (b #2a((1 2 3)
		 (4 5 6)
		 (7 8 9))))
      (with-matrixes
	  (map + A B)
	:field fixnum))
  #2A ((2 3 4) (5 6 7) (8 9 10)))

(rt:deftest (linear-combination-mismatch mat)
  (handler-case
      (with-matrixes
	  (+ (* 1 #2A ((1 2 3)))
	     (* 1 #2a ((1 2))))
	  :field fixnum)
    (t (e) (type-of e)))
  matrix-error)

(rt:deftest normalize
    (multiple-value-bind (x a)
	(regression::normalize #2A((1 0.5 1) (1 -0.5 3)) )
      (values
       (with-matrixes
	   (map round X) :field fixnum)
       (with-matrixes
	   (map round a) :field fixnum)))
  #2A ((1 1 -1) (1 -1 1))
  #2a ((1 -0 -2) (0 2 0) (0 0 1)))

(rt:deftest gamma-matrixes
    (let* ((gamma
	     (vector
	      #2A((1 0 0 0) (0 1 0 0) (0 0 -1 0) (0 0 0 -1))
	      #2A((0 0 0 1) (0 0 1 0)  (0 -1 0 0) (-1 0 0 0))
	      #2A((0 0 0 #C (0 -1)) (0 0 #C (0 1) 0) (0 #C (0 1) 0 0) (#C (0 -1) 0 0 0))
	      #2A((0 0 1 0) (0 0 0 -1) (-1 0 0 0) (0 1 0 0))))
	   (res (make-array '(4 4))))
      (dotimes (i 4)
	(dotimes (j 4)
	  (setf (aref res i j)
		(position
		 (with-matrixes (+
				 (* (aref gamma i) (aref gamma j))
				 (* (aref gamma j) (aref gamma i)))
		   :field t)
		 #(#2a((2 0 0 0) (0 2 0 0) (0 0 2 0) (0 0 0 2))
		   #2a((0 0 0 0) (0 0 0 0) (0 0 0 0) (0 0 0 0))
		   #2a((-2 0 0 0) (0 -2 0 0) (0 0 -2 0) (0 0 0 -2)))
		 :test 'equalp))))
      res)
  #2A((0 1 1 1) (1 2 1 1) (1 1 2 1) (1 1 1 2)))

(rt:deftest pauli-fixnum
    (with-matrixes (:pauli :x) :field fixnum)
  #2A((0 1) (1 0)))

(rt:deftest pauli-fixnum-type
    (type-of (with-matrixes (:pauli :x) :field fixnum))
  (simple-array fixnum (2 2)))

(rt:deftest pauli-single
    (with-matrixes (:pauli :x))
  #2A((0.0 1.0) (1.0 0.0)))

(rt:deftest pauli-single-type
    (type-of (with-matrixes (:pauli :x)))
  (simple-array single-float (2 2)))
