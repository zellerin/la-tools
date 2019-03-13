(asdf:defsystem la-tools
  :version "0"
  :description "Machine learning course exercises reworked in CL"
  :maintainer "Tomas Zellerin <zellerin@gmail.com>"
  :author "Tomas Zellerin <zellerin@gmail.com>"
  :licence "BSD-style"
  :depends-on (rt alexandria)
  :serial nil
  ;; components likely need manual reordering
  :components ((:static-file "README.org" :pathname "README.org")
	       (:file "package")
	       (:file "matrixes")
	       (:file "utils")
	       (:file "multiply-test")
	       (:file "regression")
	       (:file "regression-test")
	       (:file "matrixes-tests")
	       )
  ;; :long-description ""
  )
