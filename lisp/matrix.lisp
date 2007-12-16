(in-package :kempbasu)

(def-input-fun read-raw-matrix (file)
    "Read a file with a matrix in the dumb format;

rows cols
x11 ... x1n
... ... ...
xm1 ... xmn"
    
    (let* ((rows (read file))
	   (cols (read file))
	   (m (make-array (list rows cols))))
      (loop for i from 0 to (1- rows)
	 do (loop for j from 0 to (1- cols)
		   do (setf (aref m i j) (read file))))
      m))


(def-input-fun read-matrix-as-lists (file &optional (delimiter #\Tab))
  "Read a file as a list of lists"
  (collect-lines file #'(lambda (x)
			  (split-sequence delimiter x))))



