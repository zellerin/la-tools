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
      for (x y) = (read-from-string (concatenate 'string "(" line ")"))
      collect x into xses
      collect y into yses
      finally (return (values size xses yses)))))

  (multiple-value-bind (size xlist ylist)
      (read-comma-file))


(defun get-lr-coefficients (count size xlist ylist)
  (let ((x (make-array `(,size 2) :initial-element 0s0 :element-type 'single-float))
	(y (make-array `(,size 1) :initial-element 0s0 :element-type 'single-float))
	(A (make-random-array 2 1 1s-4)))
    (loop for row from 0
	  for xval in xlist
	  and yval in ylist
	  do
	     (setf (aref x row 0) 1s0
		   (aref x row 1) xval
		   (aref y row 0) (coerce yval 'single-float)))
    (dotimes (i count)
      (print (linear-regression-iteration y A x -1.25e-4 1s0)))
    (list (aref A 0 0) (aref A 1 0))))
