#-- Classes and functions for colored python printout
# Link : https://www.geeksforgeeks.org/print-colors-python-terminal/
# Link : https://pypi.org/project/colorama/ ; colorama is cross-platform

# //--------------------------------------------
# //--------------------------------------------

# Python program to print colored text and background
'''Colors class:reset all colors with colors.reset;
two sub classes : fg for foreground and bg for background;
use as colors.subclass.colorname, i.e. colors.fg.red or colors.bg.green
Also, the generic bold, dim, underline, reverse, strike through, and invisible work with the main class i.e. colors.bold'''
class colors:
    reset='\033[0m'
    bold='\033[01m'
    dim='\033[02m'
    ital='\033[03m'
    underline='\033[04m'
    reverse='\033[07m'
    strikethrough='\033[09m'
    invisible='\033[08m'

    class fg:
        black='\033[30m'
        red='\033[31m'
        green='\033[32m'
        orange='\033[33m'
        blue='\033[34m'
        purple='\033[35m'
        cyan='\033[36m'
        lightgrey='\033[37m'
        darkgrey='\033[90m'
        lightred='\033[91m'
        lightgreen='\033[92m'
        yellow='\033[93m'
        lightblue='\033[94m'
        pink='\033[95m'
        lightcyan='\033[96m'
    class bg:
        black='\033[40m'
        red='\033[41m'
        green='\033[42m'
        orange='\033[43m'
        blue='\033[44m'
        purple='\033[45m'
        cyan='\033[46m'
        lightgrey='\033[47m'

# //--------------------------------------------
# //--------------------------------------------

# class cTerm:
#     PURPLE = '\033[95m'
#     CYAN = '\033[96m'
#     DARKCYAN = '\033[36m'
#     BLUE = '\033[94m'
#     GREEN = '\033[92m'
#     YELLOW = '\033[93m'
#     RED = '\033[91m'
#     BOLD = '\033[1m'
#     UNDERLINE = '\033[4m'
#     END = '\033[0m'
#Example : print cTerm.GREEN+"ROOT module imported"+cTerm.END

# //--------------------------------------------
# //--------------------------------------------

# from termcolor import colored, cprint
# text = colored('Hello, World!', 'red', attrs=['reverse', 'blink'])
#Example : print(text)
# cprint('Hello, World!', 'green', 'on_red')
# print_red_on_cyan = lambda x: cprint(x, 'red', 'on_cyan')
# print_red_on_cyan('Hello, World!'

# //--------------------------------------------
# //--------------------------------------------

# def prRed(skk): print("\033[91m {}\033[00m" .format(skk))
# def prGreen(skk): print("\033[92m {}\033[00m" .format(skk))
# def prYellow(skk): print("\033[93m {}\033[00m" .format(skk))
# def prLightPurple(skk): print("\033[94m {}\033[00m" .format(skk))
# def prPurple(skk): print("\033[95m {}\033[00m" .format(skk))
# def prCyan(skk): print("\033[96m {}\033[00m" .format(skk))
# def prLightGray(skk): print("\033[97m {}\033[00m" .format(skk))
# def prBlack(skk): print("\033[98m {}\033[00m" .format(skk))
#Example : prCyan("Hello World, ")
# prYellow("It's")
# prGreen("Geeks")
# prRed("For")
# prGreen("Geeks")

# //--------------------------------------------
# //--------------------------------------------

# class bcolors:
#     HEADER = '\033[95m'
#     OKBLUE = '\033[94m'
#     OKGREEN = '\033[92m'
#     WARNING = '\033[93m'
#     FAIL = '\033[91m'
#     ENDC = '\033[0m'
#     BOLD = '\033[1m'
#     UNDERLINE = '\033[4m'

# //--------------------------------------------
# //--------------------------------------------

# from ColoredPrintout import bcolors
# print(bcolors.WARNING + "Warning: No active frommets remain. Continue?" + bcolors.ENDC)
# print(f"{bcolors.WARNING}Warning: No active frommets remain. Continue?{bcolors.ENDC}") [PYTHON3.6]

# //--------------------------------------------
# //--------------------------------------------

# from colorama import Fore, Back, Style
# print(Fore.RED + 'some red text')
# print(Back.GREEN + 'and with a green background')
# print(Style.DIM + 'and in dim text')
# print(Style.RESET_ALL)
# print('back to normal now')

# //--------------------------------------------
# //--------------------------------------------
