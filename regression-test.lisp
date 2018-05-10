(in-package regression)

(defun read-comma-file (file)
  (with-open-file (in  file)
    (loop
      with *readtable* = (copy-readtable)
	initially (set-macro-character
		   #\,
		   (lambda (s c) (declare (ignore c)) (read s)) nil)
      for line = (read-line in nil nil)
      while line
      for data = (read-from-string (concatenate 'string "(" line ")"))
      count 1 into size
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

(defun read-comma-file2 (file)
  (with-open-file (in  file)
    (loop
      with *readtable* = (copy-readtable)
	initially (set-macro-character
		   #\,
		   (lambda (s c) (declare (ignore c)) (read s)) nil)
      for line = (read-line in nil nil)
      while line
      for data = (read-from-string (concatenate 'string "(" line ")"))
      collect (butlast data) into xses
      collect (car (last data)) into yses
      count 1 into size
      finally
	 (return  (loop
		    with x = (make-array
			      `(,size (+ 2
				       (length (car xses))))
			      :initial-element 0s0 :element-type 'single-float)
		    and y = (make-array `(,size 1) :initial-element 0s0 :element-type 'single-float)
		    for row from 0
		    for xrow in xses
		    and yval in yses
		    for x0 = (car xrow)
		    do
		       (loop for xval in xrow
			     for col from 1
			     do
				(setf (aref x row col)
				      (coerce xval 'single-float)))

		       (setf (aref x row 0) 1s0
			     (aref x row (1+ (length xrow))) (/ 1s0 x0)
			     (aref y row 0) (coerce yval 'single-float))
		    finally (return (values y x)))))))
