(ns pca-example.core
  (:require [oz.core :as oz]
            [pca.core :as pca]
            [plotly-clj.core :as plot]
            [clojure.core.matrix :as matrix]
            [clojure.core.matrix.operators :as operators]
            [clojure.core.matrix.linear :as linear]
            [clojure.data.csv :as csv]
            [clojure.java.io :as io]
            [clojure.string :as string]
            [uncomplicate.commons.core :as commons]
            [uncomplicate.neanderthal.core :as neanderthal]
            [uncomplicate.neanderthal.native :as native]
            [uncomplicate.neanderthal.linalg :as linalg]
            [uncomplicate.fluokitten.core :as fluokitten]
            [clj-http.client :as http]
        )
  )

;;; COMMON

;;; Download the data from https://vincentarelbundock.github.io/Rdatasets/csv/datasets/iris.csv

(defn download-file
  [filename]
  (http/get filename))

(defonce iris-data
  (download-file "https://vincentarelbundock.github.io/Rdatasets/csv/datasets/iris.csv"))

(defn read-csv
  [file]
  (with-open [in-file (io/reader file)]
    (doall
     (csv/read-csv in-file))))

(defn format-row [row]
  (let [row (rest row)                  ; Ignore first column (index)
        floats (vec (map (fn [d i] (Float/parseFloat d)) row (range 0 4)))]
    floats))

(defn get-header
  [lines]
  (into [] (butlast (rest (first lines)))))

(defn format-data
  [raw-data]
  (map format-row (rest raw-data)))     ; Ignore the header row

;;; PCA LIB

;;; Stolen from https://gist.github.com/hswick/83f338107a9fb72082a0131e147b35ed

(defn gen-split-dataset
  [compressed-dataset labels species]
  (let [rows (map (fn [a b] (flatten [a b])) compressed-dataset labels)
        filter-fn (fn [specie] (filter #(= specie (nth % 2)) rows))]
    (map filter-fn species)))

(defn build-trace [data name]
  {:name name
   :type "scatter"
   :x (map first data)
   :y (map second data)
   :mode "markers"})

;;; OZ

(defn vecseq2mapseq
  [vs & ks]
  (map #(zipmap ks %) vs))

(defn gen-scatter-plot
  [data]
  {:width 800
   :height 800
   :data {:values data}
   :encoding {:x {:field "x"
                  :scale {:domain [4 10]}}
              :y {:field "y"
                  :scale {:domain [4 10]}}
              :color {:field "species" :type "nominal"}}
   :mark "circle"})

;;; CORE.MATRIX

(defn ds2array
  [ds]
  (matrix/array (into [] ds)))

;; See core.matrix.stats/mean
(defn emean
  [v]
  (/ (matrix/esum v) (matrix/ecount v)))

(defn colmeans
  [m]
  (into [] (map emean (matrix/slices m 1))))

(defn subtract-means
  [m]
  (let [means (colmeans m)]
    (into [] (map #(into [] (matrix/sub % means)) m))))

(defn covariance
  [m]
  (let [m-trans (matrix/transpose m)]
    (matrix/mmul m-trans m)))

(defn get-eigenvecs
  [m]
  (:Q (linear/eigen m)))

;;; NEANDERTHAL

(defn ds2dge
  [ds]
  (let [nrows (count ds)
        ncols (count (first ds))]
    (native/dge nrows ncols (apply concat ds) {:layout :row})))

(defn nmean
  [v]
  (/ (neanderthal/sum v) (neanderthal/dim v)))

(defn ncolmeans
  [m]
  (native/dv (fluokitten/fmap nmean (neanderthal/cols m))))

(defn idvec
  [n]
  (fluokitten/fmap (fn ^double [^double x] (inc x)) (native/dv n)))

(defn nsubtract-colmeans
  [m]
  (let [numrows (neanderthal/mrows m)
        id (idvec numrows)
        means (ncolmeans m)
        mmeans (neanderthal/rk id means)]
    (neanderthal/axpy -1 mmeans m)))

(defn nadd-colmeans
  [m means]
  (let [numrows (neanderthal/mrows m)
        id (idvec numrows)
        mmeans (neanderthal/rk id means)]
    (neanderthal/axpy 1 mmeans m)))

(defn ncovariance
  [m]
  (neanderthal/mm (neanderthal/trans m) m))

(defn get-eigens
  [m]
  (let [mc (neanderthal/copy m)
        numrows (neanderthal/mrows mc)
        numcols (neanderthal/ncols mc)
        eigenvals (native/dge numrows 2)
        left-eigenvecs (native/dge numrows numcols)
        right-eigenvecs (native/dge numrows numcols)]
    (linalg/ev! mc eigenvals left-eigenvecs right-eigenvecs)
    {:qr-factors mc
     :eigenvals eigenvals
     :left-eigenvecs left-eigenvecs
     :right-eigenvecs right-eigenvecs
     }))

(defn format-axis-label
  [hdr ev]
  (string/join " + " (map #(str (format "%.3f" %1) " * " %2) ev hdr)))

(defn format-axis-labels
  [hdr colevecs]
  (map #(format-axis-label hdr %) colevecs))

(defn oz-data
  [coleigenvecs xs ys lines]
  (let [header (get-header lines)
        labels (format-axis-labels header coleigenvecs)
        xlabel (first labels)
        ylabel (second labels)
        zs (rest (map #(nth % 5) lines))
        ds (map (fn [x y z] {:x x :y y :species z}) xs ys zs)]
    {:width 800
     :height 800
     :data {:values ds}
     :encoding {:x {:field "x"
                    :axis {:title xlabel}
                    }
                :y {:field "y"
                    :axis {:title ylabel}
                    }
                :color {:field "species" :type "nominal"}}
     :mark "circle"}))

(comment
  ;; COMMON
  (def iris-lines (->> iris-data
                       :body
                       string/split-lines
                       (map #(-> %
                                 (string/replace #"\"" "")
                                 (string/split #",")
                                 ))))
  (def header (get-header iris-lines))
  (def dataset (format-data iris-lines))
  (oz/start-plot-server!)
  ;; PCA LIB
  (def covar (pca/covariance dataset))
  (def iris-components (pca/principal-components covar 2))
  (def compressed-dataset (map (partial pca/compress iris-components) dataset))
  (def labels (map #(get % 5) iris-lines))
  (def species ["setosa" "versicolor" "virginica"])
  (def split-dataset (gen-split-dataset compressed-dataset labels species))
  (def iris-trace (map build-trace split-dataset species))
  (plot/plot 800 800 iris-trace {})
  ;; OZ
  (def iris-data (vecseq2mapseq (apply concat split-dataset) :x :y :species))
  (def iris-plot (gen-scatter-plot iris-data))
  (oz/v! line-plot)
  (oz/v! iris-plot)
  ;; CORE.MATRIX
  (def iris-array (ds2array dataset))
  (def meansadj-iris-array (subtract-means iris-array))
  (def cov-iris-array (covariance meansadj-iris-array))
  (def eigenvecs-array (get-eigenvecs cov-iris-array))
  (def top2evecs-array (matrix/select eigenvecs-array :all [0 1]))
  (def reduced-iris-array (matrix/mmul meansadj-iris-array top2evecs-array))
  (def oz-iris-data-array (oz-data (matrix/columns top2evecs-array) (first (matrix/columns reduced-iris-array)) (second (matrix/columns reduced-iris-array)) iris-lines))
  (oz/v! oz-iris-data-array)
  ;; NEANDERTHAL
  (def iris-dge (ds2dge dataset))
  (def iris-means-vec (ncolmeans iris-dge))
  (def meansadj-iris-dge (nsubtract-colmeans iris-dge))
  (def cov-iris-dge (ncovariance meansadj-iris-dge))
  (def eigenstuff (get-eigens cov-iris-dge)) ; Are the results of eigenstuff always sorted?
  (def eigenvecs (:left-eigenvecs eigenstuff))
  (def top2evecs (neanderthal/submatrix eigenvecs 4 2)) ; Note this does not match docs
  (def reduced-iris-dge (neanderthal/mm meansadj-iris-dge top2evecs))
  ;(def iris-meansunadj-vec (native/dv (map #(neanderthal/dot iris-means-vec %) (neanderthal/cols top2evecs))))
  ;(def reduced-meansunadj-iris-dge (nadd-colmeans reduced-iris-dge iris-meansunadj-vec))
  (def oz-iris-data-dge (oz-data (neanderthal/cols top2evecs) (neanderthal/col reduced-iris-dge 0) (neanderthal/col reduced-iris-dge 1) iris-lines))
  (oz/v! oz-iris-data-dge)
  ;; TEST
  (def toym (native/dge 5 2 [2.5 3 3.25 3.75 5.1 4.9 5.8 5.7 8.1 7.2] {:layout :row}))
  (def toymeans (ncolmeans m))
  (def toymeansadj (nsubtract-colmeans m))
  (def toycovar (ncovariance tmeansadj))
  (def toyeigens (get-eigens tcovar))
  (def toyeigenvecs (:left-eigenvecs teigens))
  (def toyprojected (neanderthal/mm tmeansadj teigenvecs))
  (def toymeansunadjusted (nadd-colmeans toyprojected toymeans))
  ;; NOTES
  ;; The Neanderthal PCA seems to match https://www.math.umd.edu/~petersd/666/html/iris_pca.html
  ;; covar (PCA LIB) does not match either cov-iris-array or cov-iris-dge. The latter 2 match.
  ;; Likely has to do with the means adjusting.
  )
