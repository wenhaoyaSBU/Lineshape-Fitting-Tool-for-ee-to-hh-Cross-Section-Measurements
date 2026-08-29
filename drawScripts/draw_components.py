import os
import ROOT
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter
import matplotlib.ticker as mticker
import argparse

parser = argparse.ArgumentParser(description='Draw the fit result with components and pulls.')
parser.add_argument("--ver", required = True, help="Specify the version of the output fig")
parser.add_argument("--it", required = True, help="Specify the iteration number of the output fig")
parser.add_argument("--method", choices=["gen", "wgt"], default="wgt",
                    help="Specify the method used: 'gen' for MC generating or 'wgt' for MC weighting.")
parser.add_argument("--addbr", action='store_true', help="Whether to include the Br pull in the plot.", default=False)

parser.add_argument("--dof", help="Manually speficy degrees of freedom", default=22, type=int)

postfix = parser.parse_args().method

version = parser.parse_args().ver
iteration = parser.parse_args().it
outDir = f"../output/figs/{version}/iteration{iteration}/"
# Make the directory if it doesn't exist
os.makedirs(outDir, exist_ok=True)

# draw_gammagg = True
addbr = parser.parse_args().addbr
dof = parser.parse_args().dof
if addbr:
    dof += 1

def qfCalc(MM, mo = 0.78266, mpi = 0.1349766):
    if MM < mo + mpi:
        return 0
    else:
        return np.sqrt((MM**2 - (mo + mpi)**2) * (MM**2 - (mo - mpi)**2)) / (2 * MM)

def BrCalc(C1, FF, MM, GG, GGee):
    return (C1**2 * FF**2 * GGee * qfCalc(MM)**3) / (MM**4 * GG)



B1 = 0.98823
# B2 = 0.893
B2 = 0.892

xlow = 3.048
xhigh = 3.122

fittedFile = "../output/fittedParameters.txt"
lines = []
params_postFit     = []
params_postFit_err = []
with open(fittedFile, 'r') as file:
    lines = file.readlines()
    for line in lines:
        if line.strip() and not line.startswith("#"):
            parts = line.split()
            params_postFit.append(float(parts[0]))
            params_postFit_err.append(float(parts[1]))
fittedChi2 = lines[-1].split()[0]
# print(f"Fitted chi2: {fittedChi2}")
# Round to 2 decimal places
fittedChi2_rounded = round(float(fittedChi2), 2)
print(f"Fitted chi2: {fittedChi2_rounded}")

fittedMfile = "../output/fittedM.txt"
MdevsAbs, Mfitted = [], []
MerrPre = []
MerrPost = []
with open(fittedMfile, 'r') as file:
    lines = file.readlines()
for line in lines:
    if line.strip():
        parts = line.split()
        MdevsAbs.append(float(parts[4]))
        Mfitted.append(float(parts[0]))
        MerrPre.append(float(parts[3]))
        MerrPost.append(float(parts[1]))
MdevsRel = [abs(MdevsAbs[i] / Mfitted[i]) for i in range(len(Mfitted))]
Mpull = [(MdevsAbs[i] / MerrPre[i] if MerrPre[i] != 0 else 0) for i in range(len(Mfitted))]
MpullErr = [(MerrPost[i] / MerrPre[i] if MerrPre[i] != 0 else 0) for i in range(len(Mfitted))]


# -----------------------------
# Load theoretical calculations
# -----------------------------
theoPoints              = np.loadtxt('../output/getpoint.txt')
theoPoints_con          = np.loadtxt('../output/getpoint_dressed_con.txt')  # continuum component
theoPoints_res          = np.loadtxt('../output/getpoint_dressed_res.txt')  # resonance component
theoPoints_int          = np.loadtxt('../output/getpoint_dressed_int.txt')  # interference component
theoPoints_res_gammagg  = np.loadtxt('../output/getpoint_dressed_res_gammagg.txt')  # resonance gammagg component

W_range             = theoPoints[:, 0]
W_range_con         = theoPoints_con[:, 0]
W_range_res         = theoPoints_res[:, 0]
W_range_int         = theoPoints_int[:, 0]
W_range_res_gammagg = theoPoints_res_gammagg[:, 0]

y_range              = theoPoints[:, 1]
y_range_con          = theoPoints_con[:, 1]
y_range_res          = theoPoints_res[:, 1]
y_range_int          = theoPoints_int[:, 1]
y_range_res_gammagg  = theoPoints_res_gammagg[:, 1]

# If W ranges are not the same, terminate the program
if not (np.array_equal(W_range, W_range_con) and
        np.array_equal(W_range, W_range_res) and
        np.array_equal(W_range, W_range_int)):
    raise ValueError("W ranges in theoretical calculation files do not match.")

# ------
#  Data
# ------
x3, y3, yerr3 = [], [], []

with open('../data/para.txt', 'r') as file:
    lines = file.readlines()


rel_sys_corrs = []
other_lines = lines[1:]
# part indices: 0 1    2       3    4       5   6    7       8  9     10        11
#               W BEMS BEMSerr Nobs NobsErr Eff Lumi LumiErr Se SeErr SysUncorr SysCorr
for line in other_lines:
    if line.strip():
        parts = line.split()
        parts = [float(p) for p in parts]
        # rel_sys_corr
        rel_sys_corrs.append(parts[11])
        # x3 taken from fittedM.txt below
        ydata = parts[3] / parts[5] / parts[6] / B1 / B1 / B2
        y3.append(ydata)
        err_stat = parts[4] / parts[5] / parts[6] / B1 / B1 / B2
        err_syst_unc = parts[10] / 100.0 * ydata
        err_syst_cor = parts[11] / 100.0 * ydata
        errTotal = np.sqrt(err_stat**2 + err_syst_unc**2 + err_syst_cor**2)
        yerr3.append(errTotal)
# San check: all the rel_sys_corrs should be the same, if not, terminate the program
if not all(corr == rel_sys_corrs[0] for corr in rel_sys_corrs):
    raise ValueError("Relative systematic correlation values in para.txt are not all the same.")
rel_sys_corr = rel_sys_corrs[0]

with open('../output/fittedM.txt', 'r') as file:
    lines = file.readlines()

for line in lines:
    if line.strip():
        parts = line.split()
        x3.append(float(parts[0]))

# Other nuisance parameters
nuisances_names = [r"$M$", r"$\Gamma$", r"$\Gamma_{ee}$", r"$S_E$", r"$f$"]
nuisances_nominals =   [3.096900, 92.6e-6, 5.971e-2*92.6e-6, 0.000900, 1.0]
nuisances_errors_pre = [0.000006, 1.70e-6, 1.0574368e-07,    0.000030, rel_sys_corr/100.0]
if addbr:
    nuisances_names.append(r"$\mathcal{B}$")
    nuisances_nominals.append(5.13e-4)
    nuisances_errors_pre.append(0.14e-4)
nuisances_indices = [0, 6, 7, 8, 9]  # indices in params_postFit corresponding to the nuisance parameters
# 0: M, 6: Γ, 7: Γee, 8: SE, 9: f.
nuisances_postFit =     [params_postFit[i] for i in nuisances_indices]
nuisances_errors_post = [params_postFit_err[i] for i in nuisances_indices]
if addbr:
    print(f"Adding Br pull.")
    BrResultFile = "fitted_Br.txt"
    with open(BrResultFile, 'r') as file:
        lines = file.readlines()
    line0 = lines[0].split()
    Br_postFit = float(line0[0])
    Br_postFit_err = float(line0[1])
    nuisances_postFit.append(Br_postFit)
    nuisances_errors_post.append(Br_postFit_err)
nuisances_pulls = [(nuisances_postFit[i] - nuisances_nominals[i]) / nuisances_errors_pre[i] for i in range(len(nuisances_names))]
nuisances_pulls_err = [nuisances_errors_post[i] / nuisances_errors_pre[i] for i in range(len(nuisances_names))]
# Print out the pulls and errors for test
for i in range(len(nuisances_names)):
    print(f"Nuisance: {nuisances_names[i]}, Pull: {nuisances_pulls[i]:.2f}, Pull Error: {nuisances_pulls_err[i]:.2f}")

# XS pulls
XSpull = [(y3[i] - np.interp(x3[i], W_range, y_range)) / yerr3[i] for i in range(len(x3))]
# XS pulls errors
XSpullErr = [yerr3[i] / yerr3[i] for i in range(len(x3))]  
# which is currently 1. **This might be updated when the theoretical uncert. are calculated**

# Check lenth vs. M
if len(x3) != len(Mfitted):
    raise ValueError("Length of x3 does not match length of Mfitted. Check fittedM.txt format.")

# ----------------------------
# Plot: main + interference inset
# ----------------------------
fig = plt.figure(figsize=(8, 8))

# Main axis
ax_main = fig.add_subplot(111)
# Adjust the size of main plot
ax_main.set_position([0.1, 0.45, 0.85, 0.5])  # [left, bottom, width, height]

# Total fit, continuum, and resonance on a log-y axis
ax_main.plot(W_range,     y_range,                     'r-',  label=r'$\sigma^{\mathrm{obs}}$ Fitting result')
ax_main.plot(W_range_con, y_range_con,                 'b--', label=r'$\sigma^{\mathrm{dressed}}$ Continuum')
ax_main.plot(W_range_res, y_range_res,                 'orange', linestyle='--', label=r'$\sigma^{\mathrm{dressed}}$ Resonance')
ax_main.plot(W_range_res_gammagg, y_range_res_gammagg, 'magenta', linestyle='--', label=r'$\sigma^{\mathrm{dressed}}$ Resonance ($\gamma^{*}gg$)')

# Data points
ax_main.errorbar(x3, y3, yerr=yerr3, fmt='.', ecolor='black', color='black',
                 label=r'$\sigma^\mathrm{mea}$ Data')

# ax_main.set_ylabel(r'$\sigma^{\mathrm{mea}}\ \mathrm{(pb)}$', fontsize=16)
ax_main.set_ylabel(r'$\sigma_{e^{+}e^{-} \to \omega\,\pi^{0}}\ \mathrm{(pb)}$', fontsize=16)
ax_main.set_xlabel(r'$W\ \mathrm{(GeV)}$', fontsize=16)
ax_main.tick_params(axis='both', labelsize=16)
ax_main.minorticks_on()
ax_main.xaxis.set_minor_locator(mticker.AutoMinorLocator(5))  # 5 subdivisions per major interval

# Use log scale for main panel
ax_main.set_yscale('log')

# Set the y limit
ax_main.set_ylim(20, 7*10**4)

# Optional limits – adjust if you like
# ax_main.set_ylim((0.1, 1e5))
ax_main.set_xlim((xlow, xhigh))

# Text for chi2/ndf (position may need adjustment)
# ax_main.text(3.111, 4e4, r'$\chi^2/\mathrm{ndf}=24.39/20$', fontsize=15, ha='center')   # MC generating
# ax_main.text(3.111, 4e4, r'$\chi^2/\mathrm{ndf}=23.97/20$', fontsize=15, ha='center')   # MC weighting
ax_main.text(3.111, 4e4, rf'$\chi^2/\mathrm{{ndf}}={fittedChi2_rounded}/{dof}$', fontsize=15, ha='center')   # Fitting result

ax_main.legend(fontsize=12, loc='upper left', frameon=False)

# Inset axis for interference term (linear scale, can go negative)
ax_inset = fig.add_axes([0.17, 0.56, 0.4, 0.21])

W_inset_min = 3.09
W_inset_max = 3.102
mask_inset_int   = (W_range_int >= W_inset_min) & (W_range_int <= W_inset_max)
mask_inset_total = (W_range     >= W_inset_min) & (W_range     <= W_inset_max)

# --- new: total (fitting result) in inset ---
ax_inset.errorbar(x3, y3, yerr=yerr3, fmt='.', ecolor='black', color='black')
ax_inset.plot(W_range[mask_inset_total],
              y_range[mask_inset_total],
              'r-', label=r'$\sigma^{\mathrm{obs}}$ Fitting result')
# also the continuum
# ax_inset.plot(W_range[mask_inset_total],
#               y_range_con[mask_inset_total],
#               'b--', label='Continuum')

# interference as before
ax_inset.plot(W_range_int[mask_inset_int],
              y_range_int[mask_inset_int],
              linestyle='-.', color='green', label=r'$\sigma^{\mathrm{dressed}}$ Interference')

ax_inset.set_xlim(W_inset_min, W_inset_max)
ax_inset.tick_params(axis='both', labelsize=10)

# ax_inset.legend(fontsize=10, loc='best', frameon=False)
ax_inset.legend(fontsize=8, loc='upper left', frameon=False)


# NEW: add subplots below main for deviations (pulls) of XS and CMS energies, and to the right of the main plot for pulls of nuisance parameters.
ax_dev_XS = fig.add_axes([0.1, 0.22, 0.85, 0.15]) # [left, bottom, width, height]
ax_dev_XS.xaxis.set_ticks_position('top')
ax_dev_XS.xaxis.set_major_formatter(NullFormatter())
ax_dev_XS.set_xlim((xlow, xhigh))
ax_dev_XS.plot([xlow, xhigh], [0, 0], 'r-', linewidth=1)  # horizontal line at y=0
ax_dev_XS.errorbar(x3, XSpull, yerr=XSpullErr, fmt='.', ecolor='black', color='black',
                 label='Cross Section Residual')  # pull values for cross section
ax_dev_XS.legend(fontsize=12, loc='lower left', frameon=False)
ax_dev_XS.set_ylabel(r'$\frac{f \cdot \sigma^\mathrm{mea} - \sigma^\mathrm{obs}}{\Delta \sigma^\mathrm{mea}}$', fontsize=16)
ax_dev_XS.minorticks_on()
ax_dev_XS.xaxis.set_minor_locator(mticker.AutoMinorLocator(5))  # 5 subdivisions per major interval

ax_dev = fig.add_axes([0.1, 0.05, 0.85, 0.15]) # [left, bottom, width, height]
# Set the horizontal tickmark of the subplot to be in the upper side and hide the values at the tickmarks
ax_dev.xaxis.set_ticks_position('top')
ax_dev.xaxis.set_major_formatter(NullFormatter())
ax_dev.yaxis.set_ticks_position('left')
ax_dev.set_xlim((xlow, xhigh))
ax_dev.set_ylabel(r'$\frac{W - W^\mathrm{Prop.}}{\Delta W}$', fontsize=16)

# horizontal line at y=0
ax_dev.plot([xlow, xhigh], [0, 0], color='black', linestyle='--', linewidth=1)  
# horizontal shaded area for 1σ, 2σ
ax_dev.fill_between([xlow, xhigh], 1, -1, color='grey', alpha=0.2)
ax_dev.fill_between([xlow, xhigh], 2, -2, color='grey', alpha=0.15)


ax_dev.errorbar(x3, Mpull, yerr=MpullErr, fmt='.', ecolor='blue', color='blue',
                 label='C.O.M Energy Pull')  # pull values for CMS energy
ax_dev.legend(fontsize=12, loc='upper left', frameon=False)
ax_dev.minorticks_on()
ax_dev.xaxis.set_minor_locator(mticker.AutoMinorLocator(5))  # 5 subdivisions per major interval


# Inset axis for nuisance parameters pulls
ax_nuisance = fig.add_axes([0.78, 0.64, 0.15, 0.2])  # [left, bottom, width, height]
# Draw in a vertical pattern
# Vertical line at y=0
ax_nuisance.axvline(0, color='black', linestyle='--', linewidth=1)
ax_nuisance.errorbar(nuisances_pulls, range(len(nuisances_names)), xerr=nuisances_pulls_err,
                     fmt='o', color='dodgerblue', ecolor='dodgerblue', capsize=5, label='Nuisance Pull')
# vertical shaded area for 1σ, 2σ
ax_nuisance.fill_betweenx(range(-1, len(nuisances_names)+1), 1, -1, color='grey', alpha=0.2)
ax_nuisance.fill_betweenx(range(-1, len(nuisances_names)+1), 2, -2, color='grey', alpha=0.15)
ax_nuisance.set_ylim(-0.5, len(nuisances_names)-0.5)  # set y-limits to fit all nuisance parameters

# Change y-ticks to show nuisance parameter names
ax_nuisance.set_yticks(range(len(nuisances_names)))
ax_nuisance.set_yticklabels(nuisances_names, fontsize=12, fontweight='bold')
ax_nuisance.legend(fontsize=12, loc=[-0.35, 1], frameon=False)
ax_nuisance.set_xlabel(r'$\frac{\theta - \overline{\theta}}{\Delta \theta}$', fontsize=16, labelpad=0)


# Save figure
fig.savefig(outDir + 'fitsta_with_pull_LC_components_' + postfix + '.pdf')
fig.savefig(outDir + 'fitsta_with_pull_components_' + postfix + '.pdf')

print("Fitting result saved to " + outDir + "fitsta_with_pull_components_" + postfix + ".pdf")

# ---------------------------------------
# Save the observed cross section to file
# ---------------------------------------
with open('observed_cross_section.txt', 'w') as f:
    f.write('W(GeV) sigma_mea(pb) sigma_err(pb)\n')
    for i in range(len(x3)):
        f.write(f'{x3[i]} {y3[i]} {yerr3[i]}\n')
