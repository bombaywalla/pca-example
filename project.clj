(defproject pca-example "0.1.0-SNAPSHOT"
  :description "A small project to learn and demo several libs using the excuse of PCA."
  :url "http://example.com/FIXME"
  :license {:name "MIT License"
            :url "https://opensource.org/licenses/MIT"}
  :dependencies [[org.clojure/clojure "1.9.0"]
                 [metasoarous/oz "1.3.1"]
                 [hswick/pca "0.1.0"]
                 [hswick/plotly-clj "0.1.3"]
                 [org.clojure/data.csv "0.1.4"]
                 [net.mikera/core.matrix "0.62.0"]
                 [net.mikera/vectorz-clj "0.47.0"]
                 [uncomplicate/commons "0.4.0"]
                 [uncomplicate/neanderthal "0.18.0"
                  ;; The exclusions are there since
                  ;; there is no MacOS release for these so far
                  :exclusions [org.jcuda/jcuda-natives
                               org.jcuda/jcublas-natives]]
                 [uncomplicate/fluokitten "0.6.1"]
                 [clj-http "3.8.0"]
                 ]
  :main pca-example.core
  )
