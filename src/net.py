import pandas as pd
import numpy as np
import torch
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt


part = 413
# part2pops = {
#     102: ('ESN', 'GBR', 'JPT', 'LWK', 'PEL'),
#     207: ('CEU', 'ESN', 'GBR', 'MSL', 'PEL'),
#     312: ('FIN', 'GWD', 'LWK', 'MXL', 'PEL'),
#     413: ('GBR', 'KHV', 'MSL', 'PEL', 'YRI'),
# }
# POPS = part2pops[part]

POPS = ('GBR', 'KHV', 'MSL')
pop_num = {pop: num for num, pop in enumerate(POPS)}
num_pop = {num: pop for num, pop in enumerate(POPS)}
n_pop = len(POPS)
n = int(n_pop * (n_pop - 1) / 2)

df = pd.read_csv(f'sub_features_{part}.csv', index_col=0)
df['TARGET'] = df['REAL'].apply(lambda x: pop_num[x])
df = df.drop('REAL', axis=1)
print(df.head())

X_train, X_test, y_train, y_test = train_test_split(
    df.values[:, :n],
    df['TARGET'].values,
    test_size=0.3,
    shuffle=True
)
#
X_train = torch.FloatTensor(X_train)
X_test = torch.FloatTensor(X_test)
y_train = torch.LongTensor(y_train)
y_test = torch.LongTensor(y_test)


class AncestryNet(torch.nn.Module):

    def __init__(self, n_hidden_neurons):
        super().__init__()

        self.fc1 = torch.nn.Linear(n, n_hidden_neurons)
        self.activ1 = torch.nn.Sigmoid()
        self.fc2 = torch.nn.Linear(n_hidden_neurons, n_hidden_neurons)
        self.activ2 = torch.nn.Sigmoid()
        self.fc3 = torch.nn.Linear(n_hidden_neurons, n_pop)
        self.sm = torch.nn.Softmax(dim=1)

    def forward(self, x):
        x = self.fc1(x)
        x = self.activ1(x)
        x = self.fc2(x)
        x = self.activ2(x)
        x = self.fc3(x)
        return x

    def inference(self, x):
        x = self.forward(x)
        x = self.sm(x)
        return x


ancestry_net = AncestryNet(7)

loss = torch.nn.CrossEntropyLoss()

optimizer = torch.optim.Adam(ancestry_net.parameters(), lr=3*1.0e-3)

batch_size = 10

for epoch in range(5000):
    order = np.random.permutation(len(X_train))
    for start_index in range(0, len(X_train), batch_size):
        optimizer.zero_grad()

        batch_indexes = order[start_index:start_index + batch_size]

        x_batch = X_train[batch_indexes]
        y_batch = y_train[batch_indexes]

        preds = ancestry_net.forward(x_batch)

        loss_value = loss(preds, y_batch)
        loss_value.backward()

        optimizer.step()

    if epoch % 100 == 0:
        test_preds = ancestry_net.forward(X_test)
        test_res = test_preds.argmax(dim=1)
        precision = (test_res == y_test).float().mean()
        print(precision)
        if precision > 0.83:
            results = pd.DataFrame(test_preds.detach().numpy())
            results['real'] = y_test
            results['pred'] = test_res
            results.to_csv(f'results_part{part}.csv', index=False)
            print(test_preds)
            print(y_test)
            break


plt.rcParams['figure.figsize'] = (12, 10)

n_classes = 3
plot_colors = ['g', 'orange', 'black']
plot_step = 0.02

x_min, x_max = X_train[:, 0].min() - 1, X_train[:, 0].max() + 1
y_min, y_max = X_train[:, 1].min() - 1, X_train[:, 1].max() + 1
z_min, z_max = X_train[:, 2].min() - 1, X_train[:, 2].max() + 1

xx, yy = torch.meshgrid(torch.arange(x_min, x_max, plot_step),
                        torch.arange(y_min, y_max, plot_step))

preds_class = preds.data.numpy().argmax(axis=1)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

n = 100

for i, color in zip(range(n_classes), plot_colors):
    indexes = np.where(y_train == i)
    ax.scatter(X_train[indexes, 0],
               X_train[indexes, 1],
               X_train[indexes, 2],
               c=color,
               label=POPS[i],
               cmap='Accent')

    ax.set_xlabel('GBR_KHV')
    ax.set_ylabel('GBR_MSL')
    ax.set_zlabel('KHV_MSL')
    ax.legend()


plt.show()



