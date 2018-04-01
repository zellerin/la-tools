
(defun test ()
  (dribble (merge-pathnames "log.txt"
			    (or *compile-file-pathname*
				*load-pathname*
				".")))
  (let ((generic-500-500 (make-array '(500 500) :initial-element 1))
	(single-500-500 (make-array '(500 500) :element-type 'single-float
					       :initial-element 1s0))
	(single-500-499 (make-array '(500 499) :element-type 'single-float
					       :initial-element 1s0)))
    (print "generic 500x500 x 500x500")
    (time (times generic-500-500 generic-500-500))
    (print "single float 500x500 x 500x500")
    (time (times single-500-500 single-500-500))
    (print "single float 500x500 x 500x499")
    (time (times single-500-500 single-500-499)))
  (dribble))
