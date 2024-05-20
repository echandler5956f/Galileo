from kan import *
# create a KAN: 2D inputs, 1D output, and 5 hidden neurons. cubic spline (k=3), 5 grid intervals (grid=5).
model = KAN(width=[2,5,1], grid=5, k=3, seed=0)

# create dataset f(x,y) = exp(sin(pi*x)+y^2)
f = lambda x: torch.exp(torch.sin(torch.pi*x[:,[0]]) + x[:,[1]]**2)
dataset = create_dataset(f, n_var=2)
dataset['train_input'].shape, dataset['train_label'].shape

# train the model
model.train(dataset, opt="LBFGS", steps=20, lamb=0.01, lamb_entropy=10.);

model = model.prune()
model(dataset['train_input'])
model.plot()

# sin appears at the top of the suggestion list, which is good!
model.suggest_symbolic(0,0,0)

# x^2 appears in the suggestion list (usually not top 1), but it is fine!
model.suggest_symbolic(0,1,0)

# exp not even appears in the list (but note how high correlation of all these functions), which is sad!
model.suggest_symbolic(1,0,0)

# let's try suggesting more by changing topk. Exp should appear in the list
# But it's very unclear why should we prefer exp over others. All of them have quite high correlation with the learned spline.
model.suggest_symbolic(1,0,0,topk=15)

model.train(dataset, opt="LBFGS", steps=20);
model.plot()

# sin appears at the top of the suggestion list, which is good!
model.suggest_symbolic(0,0,0)

# x^2 appears at the top of the suggestion list, which is good!
# But note how competitive cosh and gaussian are. They are also locally quadratic.
model.suggest_symbolic(0,1,0)

# exp appears at the top of the suggestion list, which is good!
model.suggest_symbolic(1,0,0)

# now let's replace every activation function with its top 1 symbolic suggestion. This is implmented in auto_symbolic()
model.auto_symbolic()

# if the user wants to constrain the symbolic space, they can pass in their symbolic libarary
# lib = ['sin', 'x^2', 'exp']
# model.auto_symbolic(lib=lib)

model.train(dataset, opt="LBFGS", steps=20);
model.plot()

# obtaining symbolic formula
formula, variables = model.symbolic_formula()
print(formula[0])