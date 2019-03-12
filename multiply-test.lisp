(in-package linear-algebra)

(defmacro time-in-ms-as-real (&body body)
  `(let ((start (get-internal-run-time)))
     ,@body
     (/ (- (get-internal-run-time) start) 0.001 internal-time-units-per-second)))

(deftest direct
    (with-matrixes
	(*
	 #2A((2 1 4 5)
	     (0 -3 -1 7)
	     (6 2 9 -8))

	 #2a((3  6)
	     (-1 1)
	     (5 0)
	     (2 -4))))

  #2A((35.0 -7.0)
      (12.0 -31.0)
      (45.0 70.0)))

(deftest transposed
  (with-matrixes

      (* (transpose
	  #2A((2 0 6)
	      (1 -3 2)
	      (4 -1 9)
	      (5 7 -8)))
	 #2a((3  6)
	     (-1 1)
	     (5 0)
	     (2 -4)))
    :field fixnum)

  #2A((35 -7)
	 (12 -31)
	 (45 70)))

(deftest (linear-combination 1)
    (with-matrixes
	(+ (* 3 #.(make-array '(3 3) :initial-element 1))
	   (* 2 #2a((1 2 3)
		   (4 5 6)
		   (7 8 9))))
      :field fixnum)
  #2A ((5 7 9) (11 13 15) (17 19 21)))

(deftest (linear-combination-mismatch)
    (handler-case
	(with-matrixes
	    (+  #2A ((1 2 3)) #2a ((1 2))))
      (t (e) (type-of e)))
  linear-algebra::matrix-error)
