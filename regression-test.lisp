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
      finally
	 (return  (loop
		    with x = (make-array
			      `(,size (1+
				       (length (car xses))))
			      :initial-element 0s0 :element-type 'single-float)
		    and y = (make-array `(,size 1) :initial-element 0s0 :element-type 'single-float)
		    for row from 0
		    for xrow in xses
		    and yval in yses
		    do
		       (loop for xval in xrow
			     for col from 1
			     do
				(setf (aref x row col)
				      (coerce xval 'single-float)))

		       (setf (aref x row 0) 1s0
			     (aref y row 0) (coerce yval 'single-float))
		    finally (return (values y x)))))))
