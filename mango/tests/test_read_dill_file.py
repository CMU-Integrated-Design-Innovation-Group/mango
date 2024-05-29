import dill

with open('./test_outputs/shape_annealing_test.aj1', 'rb') as f:
    data = dill.load(f)

print(data)