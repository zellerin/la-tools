(in-package regression)

(defun read-comma-file (file)
  (with-open-file (in  file)
    (loop
      with *readtable* = (copy-readtable)
	initially (set-macro-character
		   #\,
		   (lambda (s c) (declare (ignore c)) (read s)) nil)
      for line = (read-line in nil nil)
      for size from 0
      while line
      for data = (read-from-string (concatenate 'string "(" line ")"))
      collect (butlast data) into xses
      collect (car (last data)) into yses
      finally (return (values size xses yses)))))

(defun get-lr-coefficients (count size xlist ylist &key (sigma -1.25e-4) (rho 1s0))
  (let ((x (make-array `(,size (1+
				  (length (car xlist)))) :initial-element 0s0 :element-type 'single-float))
	(y (make-array `(,size 1) :initial-element 0s0 :element-type 'single-float))
	(A (make-random-array (1+ (length (car xlist))) 1 1s-4)))
    (loop for row from 0
	  for xrow in xlist
	  and yval in ylist
	  do
	     (loop for xval in xrow
			    for col from 1
			    do
			       (setf (aref x row col)
				     (coerce xval 'single-float)))

	     (setf (aref x row 0) 1s0
		   (aref y row 0) (coerce yval 'single-float)))
    (dotimes (i count)
      (print (linear-regression-iteration y A x sigma rho)))
    (make-array (array-dimension A 0)
		:element-type 'single-float
		:displaced-to A)))
