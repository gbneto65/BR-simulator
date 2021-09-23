"""
Created on Fri Sep 17 13:35:22 2021
@author: BRGUBO
"""

# Bacillus Recovery Simulator
# Sept. 2021
#

import numpy as np
import sys
import scipy.stats as st
import matplotlib.pyplot as plt
import pyfiglet

plt.style.use('seaborn-bright')
#plt.style.use('bmh')


# std dev. calculation
def calc_std_dev(mean, var):
    std_dev = mean * var
    return np.median(std_dev)

# generate list with random normal numbers
def random_normal (mean,sd,n):
    random = np.random.normal(mean, sd, n)
    return random

def calc_cv (mean, std):
    std= mean - mean* (1-std)
    cv = std / mean * 100
    return np.median(cv)

def intro():
    print('\n\n***************************************************************\n')
    text_title = pyfiglet.figlet_format("Bacillus Recovery Simulator")
    print(text_title)
    print('\n*****************************************************************')

    print('Python ' + sys.version)
    print('\n*****************************************************************')
    print('User Inputs :\n')
    print(f'Probiotic CFU / g: {format(probiotic_conc_g, "10.2E")}')
    print(f'Probiotic overblending: {overblending_perc} % -- CV: {overblending_var * 100} %')
    print(f'Inclusion rate (g/ton): {inclusion_g_ton} g-- CV: {inclusion_var * 100} %')
    print(f'Mixing error  (CV): {blending_var*100} %')
    print(f'CFU losses after pelletzing process (%): {pellet_losses_perc * 100} - CV: {pellet_losses_var * 100} %')
    print(f'Analitical Methodology variability (%): {analitic_method_var*100}\n')
    print(f'Number of runs in the simulation model: {n_repetition}\n')
    print('\n*****************************************************************')
    
#############################################################################
# input variables

probiotic_conc_g = 3.2E9
inclusion_g_ton = 500

inclusion_var = .1 # CV on inclusion


blending_var = .1  # CV on blending (mixing)
blending_losses = 0 # ZERO for NO LOSSES - percentage CFU reduction during the pelletizing process

pellet_losses_perc = .3 # percentage CFU reduction during the pelletizing process
pellet_losses_var = .1 # CV on losses after pelletizing
analitic_method_var = .2 # CV in the analitical method 

acceptance_range = .3 # actual defined the variation of the acceptance levels (actually +- 30% or .3)

overblending_perc = .3 # rate of cfu overblending in the product
overblending_var = .1 # CV accepted by QC on overblending in final product


n_repetition = 100000 # number of repetition

#############################################################################
intro() # print the introduction screen



label_cfu_gram_feed = probiotic_conc_g * inclusion_g_ton / 1E6


lower_accepted_value = label_cfu_gram_feed * (1-acceptance_range)
higher_accepted_value = label_cfu_gram_feed * (1+ acceptance_range)
print('\nRange of accepted CFU levels\n')
print(f'Lower = {format(lower_accepted_value, "10.2E")}')
print(f'Average = {format(label_cfu_gram_feed, "10.2E")}')
print(f'Higher = {format(higher_accepted_value, "10.2E")}')


#z_score = z_calc(normal_curve_area)
# product final CFU considering overblending

probiotic_mean_cfu_after_overblending =  probiotic_conc_g * (1 + overblending_perc) # theorical CFU after blending


lower_factor = (1 + overblending_perc) * (1 - overblending_var)
higher_factor = (1 + overblending_perc) * (1 + overblending_var)


lower= probiotic_mean_cfu_after_overblending * (1 - overblending_var)
upper = probiotic_mean_cfu_after_overblending * (1 + overblending_var)
mu = probiotic_mean_cfu_after_overblending
sigma = probiotic_mean_cfu_after_overblending - lower

probiotic_final_mean_cfu = st.truncnorm.rvs(
    (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma, size=n_repetition)


#probiotic_final_mean_cfu = probiotic_conc_g * factor_inclusion_overblending
print(f'\nprobiotic_final_mean_cfu = {format(np.mean(probiotic_final_mean_cfu), "10.2E")}')
#print(f'probiotic_final_std_cfu {format(np.std(probiotic_final_mean_cfu), "10.2E")}')

cv_probiotic_final_cfu = np.std(probiotic_final_mean_cfu) / np.mean(probiotic_final_mean_cfu) * 100
print(f'CV CFU final product = {round(cv_probiotic_final_cfu,2)} %\n')
print('-----------------------------IN FEED ----------------------------------------')


# CFU in Feed

lower_incl= 1 * (1 - inclusion_var)
upper_incl = 1 * (1 + inclusion_var)
mu_incl = 1
sigma_incl = 1 - lower_incl

# inclusion rate was considered as a trucated normal distribution where variability define the lower and upper bonds

perc_inclusion = st.truncnorm.rvs(
    (lower_incl - mu_incl) / sigma_incl, (upper_incl - mu_incl) / sigma_incl, loc=mu_incl, scale=sigma_incl, size=n_repetition)

#plt.hist(perc_inclusion)

inclusion_g_ton_real = inclusion_g_ton * perc_inclusion

plt.hist(inclusion_g_ton_real)

CFU_feed_after_inclusion_real = probiotic_final_mean_cfu * inclusion_g_ton_real / 1E6 # calculate the CFU in feed


#print(CFU_feed_after_inclusion_real)
print(f'\nCFU after inclusion ={format(np.average(CFU_feed_after_inclusion_real),"10.2E")}')

cv_after_inclusion = np.std(CFU_feed_after_inclusion_real) / np.mean(CFU_feed_after_inclusion_real) * 100
print(f'CV after inclusion = {round(cv_after_inclusion,2)} %')



# CFU after blending
if blending_losses > 1 :
    print("ERROR - Blending Losses cannot be more than 1 (100%)")
    sys.exit()
elif blending_losses == 0:
    # no losses during blending process
    CFU_feed_after_blending_std = calc_std_dev(CFU_feed_after_inclusion_real, blending_var)
    CFU_after_blending_real = random_normal(CFU_feed_after_inclusion_real, CFU_feed_after_blending_std, n_repetition)
elif blending_losses > 0 :
    # losses during blending process
    percent_losses_blending_real = random_normal((1-blending_losses),blending_var,n_repetition)
    CFU_after_blending_real = CFU_feed_after_inclusion_real * percent_losses_blending_real
else:
    print("ERROR - Blending Losses cannot be more than 1 or less than 0 (0 - 100%)")
    sys.exit()

plt.hist(CFU_after_blending_real)






#print(CFU_after_blending_real)
print(f'\nMean CFU after blending ={format(np.average(CFU_after_blending_real), "10.2E")}')
cv_after_blending = np.std(CFU_after_blending_real) / np.mean(CFU_after_blending_real) * 100
print(f'CV after blending = {round(cv_after_blending,2)} %')

# losses pellet process
percent_losses_real = random_normal((1-pellet_losses_perc),pellet_losses_var,n_repetition)
CFU_feed_after_pellet = CFU_after_blending_real * percent_losses_real

#print(CFU_feed_after_pellet)
print(f'\nMean CFU after pellet ={format(np.average(CFU_feed_after_pellet), "10.2E")}')
cv_after_pellet = np.std(CFU_feed_after_pellet) / np.mean(CFU_feed_after_pellet) * 100
print(f'CV after pellet = {round(cv_after_pellet,2)} %')

# analitical methodology variation before pellet
CFU_feed_after_method_before_pellet_std = calc_std_dev(CFU_after_blending_real, analitic_method_var)
CFU_after_method_real_before_pellet = random_normal(CFU_after_blending_real, CFU_feed_after_method_before_pellet_std, n_repetition)
cv_after_method_before_pellet = np.std(CFU_after_method_real_before_pellet) / np.mean(CFU_after_method_real_before_pellet) * 100


# analitical methodology variation after pellet
CFU_feed_after_method_std = calc_std_dev(CFU_feed_after_pellet, analitic_method_var)
CFU_after_method_real = random_normal(CFU_feed_after_pellet, CFU_feed_after_method_std, n_repetition)





#print(CFU_after_method_real)
print('\n-----------------------------IN LAB --------------------------------------------')

print(f'CV before pellet & after method. = {round(cv_after_method_before_pellet,2)} %')

print(f'\nMean CFU after method = {format(np.average(CFU_after_method_real), "10.2E")}')
cv_after_pellet = np.std(CFU_after_method_real) / np.mean(CFU_after_method_real) * 100
print(f'CV after method = {round(cv_after_pellet,2)} %')


# estimate the rejected samples

# before pelletizing process 
perc_lower_range_before_pellet_after_method = np.sum(CFU_after_method_real_before_pellet < lower_accepted_value)/n_repetition * 100
print('\n-------------------outside_acceptable levels (OAL)-------------------------------\n')
print(f'Rejected samples (OAL) before pelletizing - lower bond: {round(perc_lower_range_before_pellet_after_method,2)} %')
perc_higher_range_before_pellet_after_method = np.sum(CFU_after_method_real_before_pellet > higher_accepted_value)/n_repetition * 100
print(f'Rejected samples (OAL) before pelletizing - higher bond: {round(perc_higher_range_before_pellet_after_method,2)} %')
perc_total_out_of_range_before_pellet_after_method = perc_lower_range_before_pellet_after_method + perc_higher_range_before_pellet_after_method
print(f'Rejected samples (OAL) before pelletizing - Total: {round(perc_total_out_of_range_before_pellet_after_method,2)} %\n')


# after pelletizing process
perc_lower_range_after_inclusion = np.sum(CFU_after_method_real < lower_accepted_value)/n_repetition * 100
print('\n-------------------outside_acceptable levels (OAL)-------------------------------\n')
print(f'Rejected samples (OAL) after pelletizing - lower bond: {round(perc_lower_range_after_inclusion,2)} %')
perc_higher_range_after_inclusion = np.sum(CFU_after_method_real > higher_accepted_value)/n_repetition * 100
print(f'Rejected samples (OAL) after pelletizing - higher bond: {round(perc_higher_range_after_inclusion,2)} %')
perc_total_out_of_range_after_inclusion = perc_lower_range_after_inclusion + perc_higher_range_after_inclusion
print(f'Rejected samples (OAL) after pelletizing - Total: {round(perc_total_out_of_range_after_inclusion,2)} %\n')




# histogram in the final product
fig, ax = plt.subplots()
ax.hist([probiotic_final_mean_cfu],
         label=['Final Product'],
         bins=30,
         alpha = .4,
         #facecolor='g',
         density = True,
         )
ax.axvline(x=probiotic_mean_cfu_after_overblending,
              color = 'r',
              linestyle = '--',
              alpha = .5,
              )
xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
y_value = (ymax - ymin)*.8
x_value = (xmax - xmin)*.9
#ax.text(lower_accepted_value * .77, y_value, f'{round(perc_lower_range_after_inclusion,1)} %')
#ax.text(x_value, y_value, f'Mean: {round(np.mean(CFU_count_after_inclusion),1)}')
#ax.text(x_value, y_value*.93, f'CV: {round(cv_after_inclusion,1)}')

ax.set_title('Histogram of predited CFU in final product',
             fontsize = 8,
             )
ax.set_xlabel('Predited CFU / g')
ax.set_ylabel('Density')
plt.savefig("product_cfu_bacillus_recovery_simul.png", dpi=120)

plt.show()


# histogram after feed inclusion
fig, ax = plt.subplots()
ax.hist([CFU_feed_after_inclusion_real,CFU_after_blending_real],
         label=['CFU after inclusion', 'CFU after blending'],
         bins=30,
         alpha = .6,
         #facecolor='r',
         density = True,
         )
ax.legend()
ax.set_xlabel('Predited CFU / g of feed')
ax.set_ylabel('Density')

ax.axvline(x=lower_accepted_value,
              color = 'r',
              linestyle = '--',
              alpha = .5,
              )
ax.axvline(x=higher_accepted_value,
              color = 'r',
              linestyle = '--',
              alpha = .5,
              )



xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
y_value = (ymax - ymin)*.8
x_value = (xmax - xmin)*.9
#ax.text(lower_accepted_value * .77, y_value, f'{round(perc_lower_range_after_inclusion,1)} %')
#ax.text(x_value, y_value, f'Mean: {round(np.mean(CFU_count_after_inclusion),1)}')
#ax.text(x_value, y_value*.93, f'CV: {round(cv_after_inclusion,1)}')

ax.set_title('Histogram of predited CFU in feed',
             fontsize = 8,
             )
plt.savefig("after_inclusion_cfu_bacillus_recovery_simul.png", dpi=120)

plt.show()


# histogram after pellet
fig, ax = plt.subplots()
ax.hist([CFU_after_blending_real,CFU_feed_after_pellet],
         label=['CFU after blending', 'CFU after pellet'],
         bins=30,
         alpha = .6,
         #facecolor='r',
         density = True,
         )
ax.legend()
ax.set_xlabel('Predited CFU / g of feed')
ax.set_ylabel('Density')
ax.axvline(x=lower_accepted_value,
              color = 'r',
              linestyle = '--',
              alpha = .5,
              )
ax.axvline(x=higher_accepted_value,
              color = 'r',
              linestyle = '--',
              alpha = .5,
              )

xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
y_value = (ymax - ymin)*.8
x_value = (xmax - xmin)*.9
#ax.text(lower_accepted_value * .77, y_value, f'{round(perc_lower_range_after_inclusion,1)} %')
#ax.text(x_value, y_value, f'Mean: {round(np.mean(CFU_count_after_inclusion),1)}')
#ax.text(x_value, y_value*.93, f'CV: {round(cv_after_inclusion,1)}')

ax.set_title('Histogram of predited CFU in feed',
             fontsize = 8,
             )
plt.savefig("after_pellet_cfu_bacillus_recovery_simul.png", dpi=120)

plt.show()

# histogram after methodoly error
fig, ax = plt.subplots()
ax.hist([CFU_feed_after_pellet, CFU_after_method_real],
         label=['CFU after pellet', 'CFU after lab method'],
         bins=30,
         alpha = .6,
         #facecolor='r',
         density = True,
         )
ax.legend()
ax.set_xlabel('Predited CFU / g of feed')
ax.set_ylabel('Density')
ax.axvline(x=lower_accepted_value,
              color = 'r',
              linestyle = '--',
              alpha = .8,
              )
ax.axvline(x=higher_accepted_value,
              color = 'r',
              linestyle = '--',
              alpha = .5,
              )

xmin,xmax = ax.get_xlim()
ymin,ymax = ax.get_ylim()
y_value = (ymax - ymin)*.8
x_value = (xmax - xmin)*.9
ax.text(lower_accepted_value * .3, y_value * .8, f'OAL = {round(perc_lower_range_after_inclusion,1)} %')

ax.set_title('Histogram of predited CFU in feed',
             fontsize = 8,
             )
plt.savefig("after_method_cfu_bacillus_recovery_simul.png", dpi=120)
plt.show()

# end of code


