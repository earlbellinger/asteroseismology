#### Deep learning for stellar structure 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

set.seed(0)

library(h2o)
localH2O <- h2o.init(nthreads=8)

train_h2o <- h2o.importFile('training_mat_r.dat') 
train_h2o <- sample(train_h2o)
test_h2o <- h2o.importFile('test_mat_r.dat')
#train_split <- h2o.splitFrame(train_h2o, ratios = 0.8)

xs <- seq(0.05, 0.35, 0.01)
xs2 <- seq(0.05, 0.35, 0.05)

prev.model <- NULL

results <- sapply(xs, function(x) {
    ii <- which(xs == x)
    model <- h2o.deeplearning(
        x = 1:ncol(test_h2o),
        y = (ncol(test_h2o)+ii),
        activation = "RectifierWithDropout",
        hidden = c(10,10,10),
        input_dropout_ratio = 0.2,
        nfolds=2, 
        epochs=3,
        training_frame=train_h2o,
        initial_weights=list(prev.model@model_id),
        initial_biases=list(prev.model@model_id))
        #training_frame = train_split[[1]],
        #validation_frame = train_split[[2]])
    print(model)
    result <- predict(model, test_h2o)
    quantiles <- data.frame(x,
                            quantile(result, 0.16), 
                            quantile(result, 0.5), 
                            quantile(result, 0.86))
    print(quantiles)
    quantiles
})

print(results)

h2o.shutdown()


library(deepnet)
training.mat <- read.table('training_mat_r.dat', header=1)
test.mat <- read.table('test_mat_r.dat', header=1)
X <- as.matrix(training.mat[,1:ncol(test.mat)])
ys <- as.matrix(training.mat[,-1:-ncol(test.mat)])
trained.nn <- nn.train(X, ys, hidden=c(10))
#trained.nn <- dbn.dnn.train(X, ys, hidden=c(10), 
#    hidden_dropout=0.05, visible_dropout=0.05, cd=1)
result <- nn.predict(trained.nn, test.mat)
apply(result, 2, mean)

