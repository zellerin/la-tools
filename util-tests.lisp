(in-package regression)

(deftest (float-sigma 0) (float-sigma 0s0) 0.5)

(deftest (float-sigma 200) (float-sigma 2s2) 1s0)

(deftest (float-sigma -200) (float-sigma -2s2)  0s0)
