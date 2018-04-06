(in-package linear-algebra)

(defmacro time-in-ms-as-real (&body body)
  `(let ((start (get-internal-run-time)))
     ,@body
     (/ (- (get-internal-run-time) start) 0.001 internal-time-units-per-second)))


(defun speed-test ()
  (dribble (merge-pathnames "log.txt"
			    (or *compile-file-pathname*
				*load-pathname*
				".")))
  (let ((generic-500-500 (make-array '(500 500) :initial-element 1))
	(single-500-500 (make-array '(500 500) :element-type 'single-float
					       :initial-element 1s0))
	(single-500-499 (make-array '(500 499) :element-type 'single-float
					       :initial-element 1s0)))
    (print "generic 500x500 x 500x500" *trace-output*)
    (time (times generic-500-500 generic-500-500))
    (print "single float 500x500 x 500x500" *trace-output*)
    (time (times single-500-500 single-500-500))
    (print "single float 500x500 x 500x499" *trace-output*)
    (time (times single-500-500 single-500-499))
    (values)))

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
