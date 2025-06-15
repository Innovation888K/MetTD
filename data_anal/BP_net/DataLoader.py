import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader, Subset
from Models import MyNet
from sklearn.model_selection import train_test_split
import numpy as np

class ExcelDataset(Dataset):
    def __init__(self, features_file=r'D:\totol_dataset.xlsx',
                 labels_file='D:\\database_detail.xlsx'):
        self.features_data = pd.read_excel(features_file)
        self.labels_data = pd.read_excel(labels_file)
        self.features = self.features_data.loc[:, '1':'111']
        self.features = self.features.T
        self.features.fillna(0, inplace=True)
        self.features = self.features.values
        self.labels = self.labels_data.iloc[:, -1].values
        self.labels, self.label_mapping = pd.factorize(self.labels)
        self.one_hot_labels = np.eye(13)[self.labels]
        assert len(self.features) == len(self.one_hot_labels), "Feature num not match"

    def __len__(self):
        return len(self.features)

    def __getitem__(self, idx):
        sample = torch.tensor(self.features[idx, :], dtype=torch.float32)
        label = torch.tensor(self.one_hot_labels[idx], dtype=torch.float32)
        return sample, label


# 使用示例
if __name__ == "__main__":
    dataset = ExcelDataset()
    train_indices = (pd.read_csv("train_indices.csv")['train_indices'] - 1).tolist()
    test_indices = (pd.read_csv("test_indices.csv")['test_indices'] - 1).tolist()
    train_dataset = Subset(dataset, train_indices)
    test_dataset = Subset(dataset, test_indices)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    input_dim = 5682
    BP = MyNet(Dem1=input_dim,
               Dem2=5 * input_dim,
               Dem3=input_dim,
               Dem4=3000,
               Dem5=2000,
               Dem6=1000,
               Dem7=13)
    train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
    num_epochs = 5  
    for epoch in range(num_epochs):
        total_loss = 0.0
        BP.train()
        for samples, labels in train_loader:
            samples = samples.to(device)
            labels = labels.to(device)
            loss = BP.backward(samples, labels)
            total_loss += loss
        print(f'Epoch [{epoch + 1}/{num_epochs}], Train Loss: {total_loss / len(train_loader):.4f}')
    BP.eval()
    test_loss = 0.0
    with torch.no_grad():
        for samples, labels in test_loader:
            samples = samples.to(device)
            labels = labels.to(device)
            predictions = BP(samples)
            loss = BP.Loss_func(predictions, labels)
            test_loss += loss.item()

    print(f'Test Loss: {test_loss / len(test_loader):.4f}')
    torch.save(BP.state_dict(), 'mynet_model.pth')