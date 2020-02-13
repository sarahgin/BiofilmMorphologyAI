"""
" The goal of this file is to load the mnist model and to use it
"""
import torch
from mnist_trial.train_mnist import Net
from torchvision import datasets, transforms


def main():
    model = Net()
    loader = torch.load("mnist_cnn.pt")
    model.load_state_dict(loader)
    print(model.eval())


    # load some test:
    test_loader = torch.utils.data.DataLoader(
        datasets.MNIST('../data', train=False, #download=True,
                       transform=transforms.Compose([
                           transforms.ToTensor(),
                           transforms.Normalize((0.1307,), (0.3081,))
                       ])),
        batch_size=1, shuffle=True)
    print(type(test_loader))

    no_cuda = False
    use_cuda = not no_cuda and torch.cuda.is_available()
    device = torch.device("cuda" if use_cuda else "cpu")
    for data, target in test_loader:
        data, target = data.to(device), target.to(device)
        print('data = {}'.format(data))
        output = model(data)
        print('output = {}'.format(output))
        break

    print('Done in main-eval')


if __name__ == '__main__':
    main()