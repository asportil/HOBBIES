/!\ This repository is a work in progress, until this banner is removed /!\

# THERMISTORS
You are a student, young engineer or just curious about thermistors and you want to use one in your hobby project.
Here I can share with you some of the knowledge I've gather by working with thermistors in different projects (research and hobbies).

However, I will focus mainly on the **accuracy** and **precision** that you can get from your acquisition system.

# What is a thermistor?

You can find very useful and comprehsive informations here : https://en.wikipedia.org/wiki/Thermistor

In simple words, the thermistor is a type of resistor whose resistance is particularly influenced by the temperature.
In general you don't want that the resitors in your circuit board to change value with temperature, however thermistors are so sensitivity to temperature variations that you can use them to reconstruct a temperature value.
The Steinhart-Hart model is a mathemetical model that help describing the relation between the temperature and the corresponding resistance value. The general form is
$$
\begin{equation}
\frac{1}{T} = A + B ln(R) + C (ln(R))^3
\end{equation}
$$
where T is the temperature, R, the thermistor resistance and A,B,C are some coefficients. 
By carefully selecting A, B, C we can get a simplified model, also known as Beta-Model.
$$
\begin{equation}
\frac{1}{T} = \frac{1}{T_0} + \frac{1}{\beta} ln(\frac{R}{R_0}) 
\end{equation}
$$

We need to focus on (at least) 2 aspects : 
-1 How do we measure the resistance of a thermistor?
-2 How do we transform it in a temperature value?
-3 How do we assess the uncertainties ?




