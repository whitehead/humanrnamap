import configparser

# Read the TOML configuration file
config = configparser.ConfigParser()
config.read('config.toml')

# Get the input variables from the configuration file if they are set, otherwise use user input
if config.has_option('general', 'AOI'):
    AOI = config.get('general', 'AOI')
else:
    AOI = input('input a title for your area of interest: ')

lst = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'M']
active = True
while active:
    if config.has_option('chromosome', 'user_chr'):
        user_chr = config.get('chromosome', 'user_chr')
    else:
        user_chr = input('input the chromosome number, X, Y, or M (hg38): ')
    if user_chr.lower() in ['x', 'y', 'm']:
        user_chr = user_chr.upper()
        active = False
    elif user_chr not in lst:
        print('invalid input.')
        active = True
    else:
        active = False

chrom = 'chr' + user_chr

# Set the maximum chromosome length based on user_chr
if user_chr == '1':
    cmax = 248956422
elif user_chr == '2':
    cmax = 242193529
# ... (rest of the chromosome length assignments)

working = True
while working:
    if config.has_option('coordinates', 'coord_opt'):
        coord_opt = config.get('coordinates', 'coord_opt')
    else:
        coord_opt = input('would you like to enter 1 or 2 sets of coordinates: ')
    if coord_opt == '1':
        running = True
        while running:
            if config.has_option('coordinates', 'start_coord1'):
                start_coord = config.get('coordinates', 'start_coord1')
            else:
                start_coord = input('input the starting coordinate (hg38): ')
            s_nums = [x for x in start_coord]
            while ',' in s_nums:
                s_nums.remove(',')
            s_str = ''.join(s_nums)
            if not s_str.isnumeric():
                print('input must consist of only numbers. try again.')
                running = True
            else:
                f_s = int(s_str)
                if cmax <= f_s:
                    print('this number is too large. try again.')
                    running = True
                else:
                    running = False
        executing = True
        while executing:
            if config.has_option('coordinates', 'end_coord1'):
                end_coord = config.get('coordinates', 'end_coord1')
            else:
                end_coord = input('input the ending coordinate (hg38): ')
            e_nums = [x for x in end_coord]
            while ',' in e_nums:
                e_nums.remove(',')
            e_str = ''.join(e_nums)
            if not e_str.isnumeric():
                print('input must consist of only numbers. try again.')
                executing = True
            else:
                f_e = int(e_str)
                if cmax < f_e:
                    print('this number is too large. try again.')
                    executing = True
                elif f_e == f_s:
                    print('this number is equal to the starting coordinate. try again.')
                    executing = True
                elif f_e < f_s:
                    print('this number is smaller than the starting coordinate. try again.')
                    executing = True
                else:
                    executing = False
                    working = False
    elif coord_opt == '2':
        # ... (rest of the code for handling two sets of coordinates)
    else:
        print("invalid input. try again.")
        working = True

managing = True
while managing:
    if config.has_option('strand', 'choice'):
        choice = config.get('strand', 'choice')
    else:
        choice = input('enter a "+" for a positive strand and a "-" for a negative strand: ')
    if choice == "+":
        strand = 'pos'
        managing = False
    elif choice == "-":
        strand = 'neg'
        managing = False
    else:
        print('invalid input. try again.')
        managing = True
