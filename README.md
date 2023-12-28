/!\ This repository is a work in progress, until this banner is removed /!\

# THERMISTORS
You are a student, young engineer or just curious about thermistors and you want to use one in your hobby project.
Here I can share with you some of the knowledge I've gather by working with thermistors in different projects (research and hobbies).

However, I will focus mainly on the **accuracy** and **precision** that you can get from your acquisition system.

We need to focus on the following aspects : 
1) What is a thermistor ?
2) How do we gather a temperature value?
3) How do we measure the resistance of a thermistor?
4) How do we assess the uncertainties ?

***What is a thermistor?***

You can find very useful and comprehsive informations here : https://en.wikipedia.org/wiki/Thermistor

In simple words, the thermistor is a type of resistor whose resistance is particularly influenced by the temperature.
In general you don't want that the resitors in your circuit board to change value with temperature, however thermistors are so sensitivity to temperature variations that you can use them to reconstruct a temperature value.

***How do we gather a temperature value?***

The Steinhart-Hart model is a mathemetical model that help describing the non-linear relationship between the temperature and the corresponding  thermistor resistance. The general form is the following

$$
\begin{equation}
  \frac{1}{T} = A + B\  ln(R) + C\   \left(ln(R) \right)^3
\end{equation}
$$

where T is the temperature, R is the thermistor resistance and A,B,C are some coefficients. 
By carefully selecting A, B, C we can get a simplified model, also known as Beta-Model.

$$
\begin{equation}
  \frac{1}{T} = \frac{1}{T_0} + \frac{1}{\beta} \  ln\left( \frac{R}{R_0} \right) 
\end{equation}
$$

where the coefficients have been selected : 

$$
\begin{equation}
\begin{aligned}
A & =   \frac{1}{T_0} + \frac{1}{\beta} ln(R) \\
B & = \frac{1}{\beta} \\
C & = 0 \\
\end{aligned}
\end{equation}
$$







