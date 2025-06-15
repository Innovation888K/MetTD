import pandas as pd
from DataLoader import ExcelDataset
import torch
from torch.utils.data import Dataset, DataLoader, Subset
from Models import MyNet
from sklearn.model_selection import train_test_split
import numpy as np


train_indices = (pd.read_csv("train_indices.csv")['train_indices']-1).tolist()
test_indices = (pd.read_csv("test_indices.csv")['test_indices']-1).tolist()

dataset = ExcelDataset()
train_dataset = Subset(dataset, train_indices)
test_dataset = Subset(dataset, test_indices)
test_loader = DataLoader(test_dataset, batch_size=32, shuffle=False)
Data_loader = DataLoader(dataset, batch_size=32, shuffle=False)
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
input_dim = 5682
BP = MyNet(Dem1=input_dim,
           Dem2=5 * input_dim,
           Dem3=input_dim,
           Dem4=3000,
           Dem5=2000,
           Dem6=1000,
           Dem7=13)
BP.load_state_dict(torch.load('mynet_model.pth'))

BP.train()
test_loss = 0.0
with torch.no_grad():
    for samples, labels in Data_loader:
        samples = samples.to(device)
        labels = labels.to(device)
        predictions = BP(samples)
        for prediction in predictions:
            print(float(prediction[0]))