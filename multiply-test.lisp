(in-package linear-algebra)

(defmacro time-in-ms-as-real (&body body)
  `(let ((start (get-internal-run-time)))
     ,@body
     (/ (- (get-internal-run-time) start) 0.001 internal-time-units-per-second)))


(deftest direct
  (equalp
   (linear-algebra:times
    #2A((2 1 4 5)
	(0 -3 -1 7)
	(6 2 9 -8))

    #2a((3  6)
	(-1 1)
	(5 0)
	(2 -4)))



   #2A((35 -7)
       (12 -31)
       (45 70)))
  t)

(deftest transposed
  (equalp

   (times-transposed
    #2A((2 0 6)
	(1 -3 2)
	(4 -1 9)
	(5 7 -8))
    #2a((3  6)
	(-1 1)
	(5 0)
	(2 -4)) )

   #2A((35 -7)
       (12 -31)
       (45 70)))
  t)
