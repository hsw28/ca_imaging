function f = f1score(ytrue, ypred)
    ytrue = double(ytrue(:));
    ypred = double(ypred(:));
    tp = sum((ytrue == 1) & (ypred == 1));
    fp = sum((ytrue == 0) & (ypred == 1));
    fn = sum((ytrue == 1) & (ypred == 0));
    prec = tp / (tp + fp + eps);
    rec  = tp / (tp + fn + eps);
    f = 2 * (prec * rec) / (prec + rec + eps);
end
