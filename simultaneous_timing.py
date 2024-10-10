import numpy as np
import matplotlib.pyplot as plt
import bilby
import psrchive
import os
import sys

from scipy.constants import c

i = sys.argv[1]
j = sys.argv[2]

i = int(i)
j = int(j)

nbins = 1024
x = np.linspace(0, 1, nbins)

noise = 0.0012


def fft_rotate(data, bins):
    """Return data rotated by 'bins' places to the left. The
       rotation is done in the Fourier domain using the Shift Theorem.
​
       Inputs:
            data: A 1-D numpy array to rotate.
            bins: The (possibly fractional) number of phase bins to rotate by.
​
       Outputs:
            rotated: The rotated data.
    """
    freqs = np.arange(data.size/2+1, dtype=np.float64)
    phasor = np.exp(complex(0.0, 2.0*np.pi) * freqs * bins / float(data.size))
    return np.fft.irfft(phasor*np.fft.rfft(data))


#Load the polarisation templates
# I don't know if these will all be the same template or not yet. For now /
# leave it like this. The templates should iterate about the j variable when ready

portrait_file = "/fred/oz002/users/mmiles/MTM_trials/J0125-2327/rm_trials/portrait_attempt/J0125-2327_wb.4p.templ_F"

port = psrchive.Archive_load(portrait_file)
port_wfreq = port.weighted_frequency(0,0,0)
port_wfreq = port_wfreq*1e6
#Attempt to remove the RM and see what changes
port.set_rotation_measure(0)

#Also turn it so it has something to fit to
port.rotate_phase(0.25)
port.remove_baseline()
port = port.get_data()
port = port[0]

## REMOVE THE /2 -- ONLY FOR TESTING
sI_template = port[0,j,:]/2
sQ_template = port[1,j,:]/2
sU_template = port[2,j,:]/2
sV_template = port[3,j,:]/2

#Rotate towards strong  - I don't think I should need this for MTM?
#weak_template = fft_rotate(weak_template,weak_template.argmax()-strong_template.argmax())


#This should load all of the data into a single archive
datas = psrchive.Archive_load('/fred/oz002/users/mmiles/MTM_trials/J0125-2327/rm_trials/data/J0125-2327_grand_T.dly_D_F_auxrm10')

rm_approx = datas.get_rotation_measure()
wfreq = datas.weighted_frequency(j,0,0)
wfreq = wfreq*1e6

datas.remove_baseline()
datas.dedisperse()
datas = datas.get_data()
#Need to keep all the data inc frequency and polarisation
datas = datas[:,:,:,:]

# def avg_profile_model(x, weak_amp, weak_phase, strong_amp, strong_phase):
 
#     weak_mode = fft_rotate(weak_amp * strong_amp * weak_template / np.max(weak_template), strong_phase + weak_phase)
#     strong_mode = fft_rotate(strong_amp * strong_template / np.max(strong_template), strong_phase)
    
#     #ratio = (np.max(weak_mode)/np.max(strong_mode))*0.15
#     #ratio = (np.max(weak_mode)/np.max(strong_mode))
#     return strong_mode + weak_mode


# def avg_profile_model(x, global_phase, sI_amp, sQ_amp, sU_amp, sV_amp):
 
#     I_mode = fft_rotate(sI_amp * sI_template / np.max(sI_template), global_phase)
#     Q_mode = fft_rotate(sQ_amp * sQ_template / np.max(sQ_template), global_phase)
#     U_mode = fft_rotate(sU_amp * sU_template / np.max(sU_template), global_phase)
#     V_mode = fft_rotate(sV_amp * sV_template / np.max(sV_template), global_phase)
    
#     #ratio = (np.max(weak_mode)/np.max(strong_mode))*0.15
#     #ratio = (np.max(weak_mode)/np.max(strong_mode))
#     return I_mode, Q_mode, U_mode, V_mode



data = datas[int(i)]
data = data[:,int(j),:]
rfft_data = np.fft.rfft(data)


rfft_sI_template = np.fft.rfft(sI_template)
rfft_sQ_template = np.fft.rfft(sQ_template)
rfft_sU_template = np.fft.rfft(sU_template)
rfft_sV_template = np.fft.rfft(sV_template)


def phasor_scale_fft(rfft_data, bins):
    """
    Add a phase gradient to the rotated FFT input
    """
    freqs = np.arange(rfft_data.size, dtype=np.float64)
    phasor = np.exp(complex(0.0, 2.0*np.pi) * freqs * bins / float(2*(rfft_data.size - 1)))
    return phasor*rfft_data

# def avg_profile_model_fdm_fast(x, weak_amp, weak_phase, strong_amp, strong_phase):
#     """
#     Model for the average profile given two modes
#     NOTE: The templates must already be defined in this script
#     """
#     weak_mode = phasor_scale_fft(weak_amp*rfft_weak_template, weak_phase)
#     strong_mode = phasor_scale_fft(strong_amp*rfft_strong_template, strong_phase)
#     return np.concatenate((np.real(strong_mode + weak_mode), np.imag(strong_mode + weak_mode)))


def avg_profile_model_fdm_fast(x, global_phase, sI_amp, sQ_amp, sU_amp, sV_amp):
    """
    Model for the average profile given two modes
    NOTE: The templates must already be defined in this script
    """


    I_mode = phasor_scale_fft(sI_amp*rfft_sI_template, global_phase)
    Q_mode = phasor_scale_fft(sQ_amp*rfft_sQ_template, global_phase)
    U_mode = phasor_scale_fft(sU_amp*rfft_sU_template, global_phase)
    V_mode = phasor_scale_fft(sV_amp*rfft_sV_template, global_phase)

    temp_comb = np.vstack([I_mode, Q_mode, U_mode, V_mode])


    #return np.concatenate((np.real(strong_mode + weak_mode), np.imag(strong_mode + weak_mode)))
    return np.concatenate((np.real(temp_comb), np.imag(temp_comb)))

def avg_profile_model_fdm_fast_RM(x, global_phase, RM, sI_amp, L_obs, sV_amp):
    """
    Model for the average profile given two modes
    NOTE: The templates must already be defined in this script
    """
    
    #L_obs = sQ_amp + 1j*sU_amp
    wavelength =  (c / wfreq) #assume freq is in Hz
    deFaradayor = np.exp(-2.0 * 1j * wavelength **2.0 * RM)
    L_deFaraday = deFaradayor * L_obs
    sQ_deF = L_deFaraday.real
    sU_deF = L_deFaraday.imag
    I_mode = phasor_scale_fft(sI_amp*rfft_sI_template, global_phase)
    Q_mode = phasor_scale_fft(sQ_deF*rfft_sQ_template, global_phase)
    U_mode = phasor_scale_fft(sU_deF*rfft_sU_template, global_phase)
    V_mode = phasor_scale_fft(sV_amp*rfft_sV_template, global_phase)

    temp_comb = np.vstack([I_mode, Q_mode, U_mode, V_mode])


    #return np.concatenate((np.real(strong_mode + weak_mode), np.imag(strong_mode + weak_mode)))
    return np.concatenate((np.real(temp_comb), np.imag(temp_comb)))

def avg_profile_model_fdm_fast_RM_QU_full(x, global_phase, RM, amp):
    """
    Model for the average profile given two modes
    NOTE: The templates must already be defined in this script
    """
    
    #L_obs = sQ_amp + 1j*sU_amp

    wavelength =  (c / wfreq) #assume freq is in Hz
    I_mode = phasor_scale_fft(amp*rfft_sI_template, global_phase)
    Q_mode = phasor_scale_fft(amp*((rfft_sQ_template*np.cos(2 * RM * wavelength ** 2)) - (rfft_sU_template*np.sin(2 * RM * wavelength ** 2))), global_phase)
    U_mode = phasor_scale_fft(amp*((rfft_sQ_template*np.sin(2 * RM * wavelength ** 2)) + (rfft_sU_template*np.cos(2 * RM * wavelength ** 2))), global_phase)
    V_mode = phasor_scale_fft(amp*rfft_sV_template, global_phase)

    temp_comb = np.vstack([I_mode, Q_mode, U_mode, V_mode])


    #return np.concatenate((np.real(strong_mode + weak_mode), np.imag(strong_mode + weak_mode)))
    return np.concatenate((np.real(temp_comb), np.imag(temp_comb)))

def avg_profile_model_fdm_fast_RM_QU_full_altfft(x, global_phase, RM, amp):
    """
    Model for the average profile given two modes
    NOTE: The templates must already be defined in this script
    """
    
    #L_obs = sQ_amp + 1j*sU_amp

    wavelength =  (c / wfreq) #assume freq is in Hz
    I_mode = phasor_scale_fft(amp*rfft_sI_template, global_phase)
    Q_mode = phasor_scale_fft(amp*np.fft.rfft((sQ_template*np.cos(2 * RM * wavelength ** 2)) - (sU_template*np.sin(2 * RM * wavelength ** 2))), global_phase)
    U_mode = phasor_scale_fft(amp*np.fft.rfft((sQ_template*np.sin(2 * RM * wavelength ** 2)) + (sU_template*np.cos(2 * RM * wavelength ** 2))), global_phase)
    V_mode = phasor_scale_fft(amp*rfft_sV_template, global_phase)

    temp_comb = np.vstack([I_mode, Q_mode, U_mode, V_mode])


    #return np.concatenate((np.real(strong_mode + weak_mode), np.imag(strong_mode + weak_mode)))
    return np.concatenate((np.real(temp_comb), np.imag(temp_comb)))

def avg_profile_model_fdm_fast_stokes_I_debug(x, global_phase, amp):
    """
    Model for the average profile given two modes
    NOTE: The templates must already be defined in this script
    """
    
    #L_obs = sQ_amp + 1j*sU_amp

    #wavelength =  (c / wfreq) #assume freq is in Hz
    I_mode = phasor_scale_fft(amp*rfft_sI_template, global_phase)
    #Q_mode = phasor_scale_fft(amp*((rfft_sQ_template*np.cos(2 * RM * wavelength ** 2)) - (rfft_sU_template*np.sin(2 * RM * wavelength ** 2))), global_phase)
    #U_mode = phasor_scale_fft(amp*((rfft_sQ_template*np.sin(2 * RM * wavelength ** 2)) + (rfft_sU_template*np.cos(2 * RM * wavelength ** 2))), global_phase)
    #V_mode = phasor_scale_fft(amp*rfft_sV_template, global_phase)

    temp_comb = I_mode
    #temp_comb = np.vstack([I_mode])


    #return np.concatenate((np.real(strong_mode + weak_mode), np.imag(strong_mode + weak_mode)))
    return np.concatenate((np.real(temp_comb), np.imag(temp_comb)))

def avg_profile_model_fdm_fast_RM_QU_full_altfft_extrappa(x, global_phase, RM, amp):
    """
    Model for the average profile given two modes
    NOTE: The templates must already be defined in this script
    """
    ## This is a really terrible model DO NOT USE
    #L_obs = sQ_amp + 1j*sU_amp

    wavelength_data =  (c / wfreq) #assume freq is in Hz
    wavelength_port = (c / port_wfreq)
    intr_ppa = 0.5*np.arctan(sU_template/sQ_template)
    psi_active = intr_ppa + (RM * ((wavelength_data**2) - (wavelength_port**2)))

    I_mode = phasor_scale_fft(amp*rfft_sI_template, global_phase)
    Q_mode = phasor_scale_fft(amp*np.fft.rfft((sQ_template*np.cos(2 * psi_active)) - (sU_template*np.sin(2 * psi_active))), global_phase)
    U_mode = phasor_scale_fft(amp*np.fft.rfft((sQ_template*np.sin(2 * psi_active)) + (sU_template*np.cos(2 * psi_active))), global_phase)
    V_mode = phasor_scale_fft(amp*rfft_sV_template, global_phase)

    temp_comb = np.vstack([I_mode, Q_mode, U_mode, V_mode])

    return np.concatenate((np.real(temp_comb), np.imag(temp_comb)))


# datarealabs = np.abs(np.real(rfft_data))
# approxmean = np.mean(datarealabs[0,-100:])
# approxstd = np.std(datarealabs[0,-100:])

data_fdm = np.concatenate((np.real(rfft_data), np.imag(rfft_data)))
noise_fdm = noise*np.sqrt(len(data_fdm)/2)

likelihood = bilby.likelihood.GaussianLikelihood(x, data_fdm,
                                                avg_profile_model_fdm_fast_RM_QU_full_altfft)

priors = dict()
#priors['weak_amp'] = bilby.core.prior.Gaussian(mu=0.116654865, sigma=0.031803366, name='weak_amp')
#priors['weak_amp'] = bilby.core.prior.Uniform(0, 1, 'weak_amp')
#priors['weak_phase'] = bilby.core.prior.Uniform(-nbins/8, nbins/8,
#                                                'weak_phase')
#priors['weak_phase'] = bilby.core.prior.Gaussian(mu=6.434, sigma=2.205031088120973, name='weak_phase')
#priors['weak_phase'] = bilby.core.prior.Uniform(-14, 0,
#                                                'weak_phase')
 

priors['RM'] = bilby.core.prior.Uniform(-100, 100, 'RM')
priors['amp'] = bilby.core.prior.Uniform(0, 10, 'amp')
#priors['sQ_amp'] = bilby.core.prior.Uniform(0, 10, 'sQ_amp')
#priors['sU_amp'] = bilby.core.prior.Uniform(0, 10, 'sU_amp')
#priors['L_obs'] = bilby.core.prior.Uniform(0, 10, 'L_obs')                       
#priors['sV_amp'] = bilby.core.prior.Uniform(0, 10, 'sV_amp')
priors['global_phase'] = bilby.core.prior.Uniform(-nbins/2, nbins/2,
                                                'global_phase')
priors['sigma'] = bilby.core.prior.Uniform(0, 10, name='sigma')

outDir = 'simultaneous_timing_{0}_Fscrunched_RMdebug_portraitRM0_altfit_altpsi'.format(i)

results = bilby.core.sampler.run_sampler(
    likelihood, priors=priors, sampler='dynesty', label='dynesty',
    nlive=100, verbose=True, resume=False, npool=1,
    outdir=outDir)

results.plot_corner()

rm = results.get_one_dimensional_median_and_error_bar('RM').median
amp = results.get_one_dimensional_median_and_error_bar('amp').median
#sqamp = results.get_one_dimensional_median_and_error_bar('sQ_amp').median
#suamp = results.get_one_dimensional_median_and_error_bar('sU_amp').median
#lobsamp = results.get_one_dimensional_median_and_error_bar('L_obs').median
#svamp = results.get_one_dimensional_median_and_error_bar('sV_amp').median
gp = results.get_one_dimensional_median_and_error_bar('global_phase').median

def fft_rotate(data, bins):
          """Return data rotated by 'bins' places to the left. The
             rotation is done in the Fourier domain using the Shift 
      Theorem.
      ​
             Inputs:
                  data: A 1-D numpy array to rotate.
                  bins: The (possibly fractional) number of phase bi
      ns to rotate by.
      ​
             Outputs:
                  rotated: The rotated data.
          """
          freqs = np.arange(data.size/2+1, dtype=np.float)
          phasor = np.exp(complex(0.0, 2.0*np.pi) * freqs * bins / float(data.size))
          return np.fft.irfft(phasor*np.fft.rfft(data))

data_si = data[0,:]
data_sq = data[1,:]
data_su = data[2,:]
data_sv = data[3,:]


# wavelength =  (c / wfreq) #assume freq is in Hz
# deFaradayor = np.exp(-2.0 * 1j * wavelength **2.0 * 0)
# L_deFaraday = deFaradayor * lobsamp
# sQ_deF = L_deFaraday.real
# sU_deF = L_deFaraday.imag


wavelength =  (c / wfreq) #assume freq is in Hz

sI_scaled_shifted = fft_rotate(amp * sI_template, gp)
sQ_scaled_shifted = fft_rotate(amp*((sQ_template*np.cos(2 * rm * wavelength ** 2)) - (sU_template*np.sin(2 * rm * wavelength ** 2))), gp)
sU_scaled_shifted = fft_rotate(amp*((sQ_template*np.sin(2 * rm * wavelength ** 2)) + (sU_template*np.cos(2 * rm * wavelength ** 2))), gp)
sV_scaled_shifted = fft_rotate(amp * sV_template, gp)


fig, axs = plt.subplots(2,2)

axs[0,0].plot(data_si, label="SI data")
axs[0,0].plot(sI_scaled_shifted, label="SI model")
axs[0,0].legend()

axs[0,1].plot(data_sq, label="SQ data")
axs[0,1].plot(sQ_scaled_shifted, label="SQ model")
axs[0,1].legend()

axs[1,0].plot(data_su, label="SU data")
axs[1,0].plot(sU_scaled_shifted, label="SU model")
axs[1,0].legend()

axs[1,1].plot(data_sv, label="SV data")
axs[1,1].plot(sV_scaled_shifted, label="SV model")
axs[1,1].legend()


fig.supylabel("Arb. Flux")
fig.supxlabel("Phase bins")


fig.savefig(outDir+'/Stokes_comparison.png'.format(i))

