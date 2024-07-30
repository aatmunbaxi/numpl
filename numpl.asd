(defsystem "cl-numpl"
  :description "Numpy-flavored incants for PetaLisp"
  :author "Aatmun Baxi"
  :depends-on
  ("petalisp"
   "alexandria")
  :components ((:file "lib/numpl")))
