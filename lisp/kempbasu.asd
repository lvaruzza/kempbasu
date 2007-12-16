(defpackage #:kempbasu-system
  (:use :cl :asdf))

(in-package #:kempbasu-system)

(defsystem :kempbasu
  :description "Statistical Differential gene expression detection system"
  :depends-on (split-sequence)
  :version "0.7"
  :author "Leonardo Varuzza <varuzza@gmail.com>"
  :license "GPLv3"
  :serial t
  :components ((:file "packages")
	       (:file "util")
	       (:file "cl-statistics")
	       (:file "IO")
	       (:file "matrix")
	       ;;(:file "gnuplot")
	       (:file "randist-gamma")
	       (:file "kemp")
	       (:file "basu")))
