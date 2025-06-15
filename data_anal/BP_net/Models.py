import pandas as pd
import torch
from torch import nn


class MyNet(nn.Module):
    def __init__(self, Dem1, Dem2, Dem3, Dem4, Dem5, Dem6, Dem7):
        super(MyNet, self).__init__()
        self.Layer1 = nn.Sequential(nn.Linear(Dem1, Dem2), nn.CELU())
        self.Layer2 = nn.Sequential(nn.Linear(Dem2, Dem3), nn.CELU())
        self.Layer3 = nn.Sequential(nn.Linear(Dem3, Dem4), nn.CELU())
        self.Layer4 = nn.Sequential(nn.Linear(Dem4, Dem5), nn.CELU())
        self.Layer5 = nn.Sequential(nn.Linear(Dem5, Dem6), nn.CELU())
        self.Layer6 = nn.Sequential(nn.Linear(Dem6, Dem7), nn.Sigmoid())
        self.Loss_func = torch.nn.MSELoss()
        self.opti = torch.optim.Adam(self.parameters(), lr=0.00001)

    def forward(self, x):
        x = self.Layer1(x)
        x = self.Layer2(x)
        x = self.Layer3(x)
        x = self.Layer4(x)
        x = self.Layer5(x)
        x = self.Layer6(x)
        if self.training:
            return x
        else:
            return (x > 0.5).float()

    def backward(self, input, Expected):
        predictions = self(input)
        loss = self.Loss_func(predictions, Expected)
        self.opti.zero_grad()
        loss.backward()
        self.opti.step()
        return loss.item()