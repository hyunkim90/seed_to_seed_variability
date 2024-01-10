###

# Set colors for plotting
phylum_colors <- c(
  "gray",'black', "#DA5724", "#5F7FC7","#508578", "#CD9BCD", "orange",
  "#5F7FC7","#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#D1A33D", "#8A7C64", "#599861","#5E738F"
)

my_color_collection <- c(
  "#CBD588", "#5F7FC7", "orange", "#AD6F3B", "#673770",
  "#D14285", "#652926", "#C84248", "#8569D5", "#5E738F",
  "#D1A33D", "#8A7C64", "#599861","#616163", "#FFCDB2",
  "#6D9F71", "#242F40",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724")

my_color_Class <- c(
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "orange", "#5F7FC7", "#CBD588", "#AD6F3B", "#673770",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',

  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724")



my_color_Order <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B",
  "#F92A82", "#ED7B84", "#5F7FC7",
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)

my_color_Family <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B",
  "#F92A82", "#ED7B84", "#5F7FC7",
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)
my_color_Family2 <- c(
  "gray",'black',"#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  '#E55812',  "#FFCDB2", "#242F40", "#6D9F71", "#CCA43B",
  "#F92A82", "#ED7B84", "#5F7FC7",
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724"
)

# up to 300
my_color_OTU <- c(
  "gray",'black', "#5F7FC7", "orange",  "#AD6F3B",
  "#673770","#D14285", "#652926", "#C84248",  "#8569D5",
  "#5E738F","#D1A33D", "#8A7C64", "#599861","#616163",
  "#FFCDB2", "#242F40", "#6D9F71",  "#CCA43B", "#F92A82",
  "#ED7B84", "#7EB77F", "#DEC4A1", "#E5D1D0", '#0E8482',
  '#C9DAEA', '#337357', '#95C623', '#E55812', '#04471C',
  '#F2D7EE', '#D3BCC0', '#A5668B', '#69306D', 'navy',
  '#1A535C', '#4ECDC4', 'orange', '#FF6B6B', "orchid1",
  'cyan2', '#FFF275', 'springgreen', '#FF3C38', '#A23E48',
  '#000000', '#CF5C36', '#EEE5E9', '#7C7C7C', '#EFC88B',

  '#2E5266', '#6E8898', '#9FB1BC', '#D3D0CB', '#E2C044',
  '#5BC0EB', '#FDE74C', '#9BC53D', '#E55934', '#FA7921',
  "#CD9BCD", "#508578", "#CBD588","#CBD588", "#5F7FC7",
  "orange",   "#AD6F3B", "#673770","#D14285", "#652926",
  "#C84248",  "#8569D5", "#5E738F","#D1A33D", "#8A7C64",
  "#599861","#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F", "#DEC4A1",
  "#E5D1D0", '#0E8482', '#C9DAEA', '#337357', '#95C623',
  '#E55812', '#04471C', '#F2D7EE', '#D3BCC0', '#A5668B',
  '#69306D', '#0E103D', '#1A535C', '#4ECDC4', '#F7FFF7',
  '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange")

my_color_otu2 <- c(
  "gray",'black', "#AD6F3B", "#673770","#D14285",
  "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D",
  "#8A7C64", "#599861", "#616163",  "#FFCDB2", "#242F40",
  "#6D9F71", "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#9FB1BC', 'springgreen', '#E2C044', '#5BC0EB', 'pink',
  "orange", "#CBD588", "#5F7FC7",
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#FFE66D', '#6699CC', '#FFF275',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#EEE5E9', '#7C7C7C', '#EFC88B', '#2E5266', '#6E8898',
  '#9FB1BC', '#D3D0CB', '#E2C044', '#5BC0EB', '#FDE74C',
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "orange")


my_color_gen <- c(
  "white",'yellow', "#AD6F3B", "#673770","#D14285",
  "#652926", "#C84248", "#8569D5", "#5E738F","#D1A33D",
  "#8A7C64", "#599861", "#616163",  "#FFCDB2", "#242F40",
  "#6D9F71", "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357',
  '#95C623', '#E55812', '#04471C', '#F2D7EE', '#D3BCC0',
  '#A5668B', '#69306D', '#0E103D', '#1A535C', '#4ECDC4',
  '#F7FFF7', '#FF6B6B', '#CD9BCD', '#6699CC', 'pink',
  '#FF8C42', '#FF3C38', '#A23E48', '#000000', '#CF5C36',
  '#9FB1BC', 'springgreen', '#E2C044', '#5BC0EB',

  "orange", "#CBD588", "#5F7FC7",
  '#9BC53D', '#E55934', '#FA7921', "#CD9BCD", "#508578", "#DA5724",
  "#CBD588", "#5F7FC7", "orange",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248",
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861",
  "#616163",  "#FFCDB2", "#242F40", "#6D9F71",
  "#CCA43B", "#F92A82", "#ED7B84", "#7EB77F",
  "#DEC4A1", "#E5D1D0", '#0E8482', '#C9DAEA', '#337357')


### Bacteria
tax_table(phy)
#### get unassigned vectors
#### get CP and MT phyloseq obj and vectors
# (1) CP
# phy.cp <- subset_taxa(phy, Class == "D_2__Chloroplast") ## just confirming
# phy.cp <- subset_taxa(phy, Order == "D_3__Chloroplast")
# vec.cp <- rownames(otu_table(phy.cp))
# length(rownames(otu_table(phy.cp))) ## 305 ASVs of CP
# vec.cp

##Latest DB
phy.cp <- subset_taxa(phy, Order == "o__Chloroplast")
vec.cp <- rownames(otu_table(phy.cp))
length(rownames(otu_table(phy.cp))) ## 62 ASVs of CP
vec.cp #41
# (2) MT
#phy.mt1 <- subset_taxa(phy, Family == "f__Mitochondria") #16 (latest DB)
phy.mt <- subset_taxa(phy, Order == "o__Rickettsiales") #17 ASVs (latest DB)
vec.mt <- rownames(otu_table(phy.mt))
tax_table(phy.mt)
length(rownames(otu_table(phy.mt))) ## 77 ASVs of MT (latest DB) 

# (3) Unassigned ans Archaea
unique(tax_table(phy)[,'Kingdom']) ## only bacteria, then no need to exclude
phy.un <- subset_taxa(phy, Kingdom %in% c("Unassigned","d__Eukaryota"))
tax_table(phy.un)
vec.un <- rownames(otu_table(phy.un))
tax_table(phy.un)
length(rownames(otu_table(phy.un))) ## 102 ASVs (Latest DB)

### exclude those vectors
phy  #864 taxa and 70 samples

### pop taxa application
## get rid of CP and MT otus
### pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### let's do it!!!

bac.clean <- pop_taxa(phy, c(vec.cp, vec.mt, vec.un))

##after clean-up
bac.clean ##644
sum(otu_table(bac.clean)) #807866
sort(colSums(otu_table(bac.clean))) #minimum 3,956 reads, maximum 41,964 reads

# checking procedure of whether the MT and CP otus are cleaned
taxa_names(subset_taxa(phy, Order == "o__Rickettsiales"))  # before 10
taxa_names(subset_taxa(phy, Order=="o__Chloroplast")) # before 305

taxa_names(subset_taxa(bac.clean , Order == "o__Rickettsiales")) # after 0
taxa_names(subset_taxa(bac.clean , Family=="f__Mitochondria")) # after 0
taxa_names(subset_taxa(bac.clean , Order == "o__Chloroplast")) # after 0



#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(phy.two)[,colnames(tax_table(phy.two))] <- gsub(tax_table(phy.two)[,colnames(tax_table(phy.two))],pattern="[A-Z]_[0-9]__",replacement="")
# phy.test4 <- phy.two %>% psmelt()
# phy.test4

tax_table(bac.clean)[,colnames(tax_table(bac.clean))] <- gsub(tax_table(bac.clean)[,colnames(tax_table(bac.clean))],pattern="[a-z]__",replacement="")

#' sample_data(phy.clean)$SampleID <- factor(sample_data(phy.clean)$SampleID, levels =target_PAB)


tax_table(bac.clean)

### filter otu with total count of 20? (in all samples)
### later we can implement
bac.clean.asv <- otu_table(bac.clean)
head(bac.clean.asv)
df.clean.asv <- data.frame(bac.clean.asv)
dim(df.clean.asv)
df.clean.asv$total <- apply(df.clean.asv, 1, sum)
head(df.clean.asv)
df.clean.asv <- tibble::rownames_to_column(df.clean.asv, var = 'ASV')


sample_names(bac.clean)

## Remove reads with over 350 bp
bac.seq <- read.fasta(file = "./Bacteria/dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")
asv_over_350bp <- attr(bac.seq[which(getLength(bac.seq)>350)], "names")
asv_less_250bp <- attr(bac.seq[which(getLength(bac.seq)<250)], "names")


bac.clean
sum(otu_table(bac.clean)) #1512840
bac.clean.ss <- pop_taxa(bac.clean,asv_over_350bp)
bac.clean.ss<- pop_taxa(bac.clean.ss,asv_less_250bp)
bac.clean.ss  ## 619 ASVs (Latest DB)
sum(otu_table(bac.clean.ss)) # 806715 (Latest DB)
sort(colSums(otu_table(bac.clean.ss))) #minimum 3,442 reads, maximum 41,964 reads

###Fungal community
#### get unassigned vectors

# (3) Unassigned
unique(tax_table(fun)[,'Kingdom'])
tax_table(fun)[,'Kingdom']
unique(tax_table(fun)[,'Kingdom'])
fun.un <- subset_taxa(fun, Kingdom %in% c("Unassigned","k__Eukaryota_kgd_Incertae_sedis","k__Viridiplantae","k__Rhizaria","k__unidentified"))
vec.un <- rownames(otu_table(fun.un))
tax_table(fun.un)
length(rownames(otu_table(fun.un))) ##  191 // 7433 ASVs (Latest DB) // 107 ASVs (allow unidentified and unassigned)

### exclude those vectors
fun  # 70 samples,  395 taxa

# sample_names(fun)
# sample_variables(fun)
# 
# ### pop taxa application
# ## get rid of CP and MT otus
# ### pop_taxa function
# pop_taxa = function(physeq, badTaxa){
#   allTaxa = taxa_names(physeq)
#   myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
#   return(prune_taxa(myTaxa, physeq))
# }
# ### Clean it up!!!
#fun.clean <- pop_taxa(fun, c(vec.un))
fun.clean <-fun
#### We will also remove the "D_3__" patterns for cleaner labels
tax_table(fun.clean)[,colnames(tax_table(fun.clean))] <- gsub(tax_table(fun.clean)[,colnames(tax_table(fun.clean))],pattern="[a-z]__",replacement="")

#' sample_data(fun.clean)$SampleID <- factor(sample_data(fun.clean)$SampleID, levels =target_PAB)

tax_table(fun.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

fun.clean    ## 395 taxa
str(fun.clean)
otu_table(fun.clean)

fun.clean  ## 3430


# ## fix it in fun.clean object!!! pop_taxa does the work
# fun.clean <- pop_taxa(fun.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(fun.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement
fun.clean.asv <- otu_table(fun.clean)
head(fun.clean.asv)
df.clean.asv <- data.frame(fun.clean.asv)
dim(df.clean.asv)
df.clean.asv$total <- apply(df.clean.asv, 1, sum)
head(df.clean.asv)
df.clean.asv <- tibble::rownames_to_column(df.clean.asv, var = 'ASV')


sample_names(fun.clean)

##### get rid of otu of less than 100 reads

fun.seq <- read.fasta(file = "./Fungi/Latest DB/dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")
asv_less_than_200bp <- attr(fun.seq[which(getLength(fun.seq)<200)], "names")
fun.clean.ss <- pop_taxa(fun.clean,asv_less_than_200bp)

fun.clean.ss ##  387
colSums(otu_table(fun.clean.ss)) ## 341856


## Designating OTU id
bac.list <- bac.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
fun.list <- fun.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))

bac.list$number <- paste0('B',1:dim(bac.list)[1])
bac.list

bac.list$OTU_id <- ifelse(is.na(bac.list$Genus),ifelse(is.na(bac.list$Family),paste0(bac.list$number,'_o_',bac.list$Order),paste0(bac.list$number,'_f_',bac.list$Family)),paste0(bac.list$number,'_',bac.list$Genus))
bac.list$OTU_id

fun.list$number <- paste0('F',1:dim(fun.list)[1])
fun.list

fun.list$OTU_id <- ifelse(is.na(fun.list$Genus),ifelse(is.na(fun.list$Family),paste0(fun.list$number,'_o_',fun.list$Order),paste0(fun.list$number,'_f_',fun.list$Family)),paste0(fun.list$number,'_',fun.list$Genus))
fun.list$OTU_id

bac.list
fun.list

otu.list <- rbind(bac.list, fun.list)
dim(otu.list)
# write.xlsx(otu.list,'otu_id_book.xlsx')


OTU_id.list <- rbind(bac.list[c('OTU','OTU_id')],fun.list[c('OTU','OTU_id')])
OTU_id.list$OTU_id




### for time series data
### Bacteria
tax_table(phy.time)
#### get unassigned vectors
#### get CP and MT phyloseq obj and vectors
# (1) CP
# phy.cp <- subset_taxa(phy, Class == "D_2__Chloroplast") ## just confirming
# phy.cp <- subset_taxa(phy, Order == "D_3__Chloroplast")
# vec.cp <- rownames(otu_table(phy.cp))
# length(rownames(otu_table(phy.cp))) ## 305 ASVs of CP
# vec.cp

##Latest DB
phy.cp <- subset_taxa(phy.time, Order == "o__Chloroplast")
vec.cp <- rownames(otu_table(phy.cp))
length(rownames(otu_table(phy.cp))) ## 20 ASVs of CP
vec.cp #20
# (2) MT
#phy.mt1 <- subset_taxa(phy, Family == "f__Mitochondria") #16 (latest DB)
phy.mt <- subset_taxa(phy.time, Order == "o__Rickettsiales") #17 ASVs (latest DB)
vec.mt <- rownames(otu_table(phy.mt))
tax_table(phy.mt)
length(rownames(otu_table(phy.mt))) ## 5 ASVs of MT (latest DB) 

# (3) Unassigned and Archaea
unique(tax_table(phy.time)[,'Kingdom']) ## only bacteria, then no need to exclude
phy.un <- subset_taxa(phy.time, Kingdom %in% c("Unassigned"))
tax_table(phy.un)
vec.un <- rownames(otu_table(phy.un))
tax_table(phy.un)
length(rownames(otu_table(phy.un))) ## 11 ASVs (Latest DB)

### exclude those vectors
phy.time  #864 taxa and 70 samples

### pop taxa application
## get rid of CP and MT otus
### pop_taxa function
pop_taxa = function(physeq, badTaxa){
  allTaxa = taxa_names(physeq)
  myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  return(prune_taxa(myTaxa, physeq))
}
### let's do it!!!

bac.time.clean <- pop_taxa(phy.time, c(vec.cp, vec.mt, vec.un))

##after clean-up
bac.time.clean ##332
sum(otu_table(bac.time.clean)) #519008
sort(colSums(otu_table(bac.time.clean))) #minimum 1876 reads, maximum 32525  reads

# checking procedure of whether the MT and CP otus are cleaned
taxa_names(subset_taxa(phy, Order == "o__Rickettsiales"))  # before 10
taxa_names(subset_taxa(phy, Order=="o__Chloroplast")) # before 305

taxa_names(subset_taxa(bac.clean , Order == "o__Rickettsiales")) # after 0
taxa_names(subset_taxa(bac.clean , Family=="f__Mitochondria")) # after 0
taxa_names(subset_taxa(bac.clean , Order == "o__Chloroplast")) # after 0



#### We will also remove the "D_3__" patterns for cleaner labels
# test success
# tax_table(phy.two)[,colnames(tax_table(phy.two))] <- gsub(tax_table(phy.two)[,colnames(tax_table(phy.two))],pattern="[A-Z]_[0-9]__",replacement="")
# phy.test4 <- phy.two %>% psmelt()
# phy.test4

tax_table(bac.time.clean)[,colnames(tax_table(bac.time.clean))] <- gsub(tax_table(bac.time.clean)[,colnames(tax_table(bac.time.clean))],pattern="[a-z]__",replacement="")

#' sample_data(phy.clean)$SampleID <- factor(sample_data(phy.clean)$SampleID, levels =target_PAB)


tax_table(bac.time.clean)

### filter otu with total count of 20? (in all samples)
### later we can implement
bac.time.clean.asv <- otu_table(bac.time.clean)
head(bac.time.clean.asv)
df.clean.asv <- data.frame(bac.time.clean.asv)
dim(df.clean.asv)
df.clean.asv$total <- apply(df.clean.asv, 1, sum)
head(df.clean.asv)
df.clean.asv <- tibble::rownames_to_column(df.clean.asv, var = 'ASV')


sample_names(bac.time.clean)

## Remove reads with over 350 bp
bac.seq <- read.fasta(file = "./Time series data/bacteria/dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")
asv_over_350bp <- attr(bac.seq[which(getLength(bac.seq)>350)], "names")
asv_less_250bp <- attr(bac.seq[which(getLength(bac.seq)<250)], "names")


bac.time.clean
sum(otu_table(bac.time.clean)) #1512840
bac.time.clean.ss <- pop_taxa(bac.time.clean,asv_over_350bp)
bac.time.clean.ss<- pop_taxa(bac.time.clean.ss,asv_less_250bp)
bac.time.clean.ss  ## 327 ASVs (Latest DB)
sum(otu_table(bac.time.clean.ss)) # 518971(Latest DB)
sort(colSums(otu_table(bac.time.clean.ss))) #minimum 1876 reads, maximum  32525  reads

###Fungal community
#### get unassigned vectors

# (3) Unassigned
unique(tax_table(fun.time)[,'Kingdom'])
tax_table(fun.time)[,'Kingdom']
unique(tax_table(fun.time)[,'Kingdom'])
fun.time.un <- subset_taxa(fun.time, Kingdom %in% c("Unassigned","k__Eukaryota_kgd_Incertae_sedis","k__Viridiplantae","k__Rhizaria","k__unidentified"))
vec.un <- rownames(otu_table(fun.time.un))
tax_table(fun.time.un)
length(rownames(otu_table(fun.time.un))) ##  191 // 7433 ASVs (Latest DB) // 107 ASVs (allow unidentified and unassigned)

### exclude those vectors
fun.time  # 70 samples,  395 taxa

# sample_names(fun.time)
# sample_variables(fun.time)
# 
# ### pop taxa application
# ## get rid of CP and MT otus
# ### pop_taxa fun.timection
# pop_taxa = fun.timection(physeq, badTaxa){
#   allTaxa = taxa_names(physeq)
#   myTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
#   return(prune_taxa(myTaxa, physeq))
# }
# ### Clean it up!!!
#fun.time.clean <- pop_taxa(fun.time, c(vec.un))
fun.time.clean <-fun.time
#### We will also remove the "D_3__" patterns for cleaner labels
tax_table(fun.time.clean)[,colnames(tax_table(fun.time.clean))] <- gsub(tax_table(fun.time.clean)[,colnames(tax_table(fun.time.clean))],pattern="[a-z]__",replacement="")

#' sample_data(fun.time.clean)$SampleID <- factor(sample_data(fun.time.clean)$SampleID, levels =target_PAB)

tax_table(fun.time.clean)
## 18. 10. 17 let's plot by otu
## showing the otu that are in the negative data otu

fun.time.clean    ## 372 taxa
str(fun.time.clean)
otu_table(fun.time.clean)

fun.time.clean  ## 372


# ## fix it in fun.time.clean object!!! pop_taxa does the work
# fun.time.clean <- pop_taxa(fun.time.clean, c('CVRG01041904.1.1229'))
# any(rownames(otu_table(fun.time.clean)) == 'CVRG01041904.1.1229') ## False!!!

### filter otu with total count of 20? (in all samples)
### later we can implement
fun.time.clean.asv <- otu_table(fun.time.clean)
head(fun.time.clean.asv)
df.clean.asv <- data.frame(fun.time.clean.asv)
dim(df.clean.asv)
df.clean.asv$total <- apply(df.clean.asv, 1, sum)
head(df.clean.asv)
df.clean.asv <- tibble::rownames_to_column(df.clean.asv, var = 'ASV')


sample_names(fun.time.clean)

##### get rid of otu of less than 100 reads

fun.time.seq <- read.fasta(file = "./Time series data/fungi/dynamic DB_time series/dna-sequences.fasta", as.string = TRUE, seqtype = "DNA")
asv_less_than_200bp <- attr(fun.time.seq[which(getLength(fun.time.seq)<200)], "names")
fun.time.clean.ss <- pop_taxa(fun.time.clean,asv_less_than_200bp)

fun.time.clean.ss ##  367
colSums(otu_table(fun.time.clean.ss)) ## 341856


## Designating OTU id
bac.time.list <- bac.time.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))
fun.time.list <- fun.time.clean.ss %>% psmelt() %>% group_by(OTU, Phylum,Class,Order,Family,Genus,Species) %>% summarise(total=sum(Abundance)) %>% arrange(desc(total))

bac.time.list$number <- paste0('B',1:dim(bac.time.list)[1])
bac.time.list

bac.time.list$OTU_id <- ifelse(is.na(bac.time.list$Genus),ifelse(is.na(bac.time.list$Family),paste0(bac.time.list$number,'_o_',bac.time.list$Order),paste0(bac.time.list$number,'_f_',bac.time.list$Family)),paste0(bac.time.list$number,'_',bac.time.list$Genus))
bac.time.list$OTU_id

fun.time.list$number <- paste0('F',1:dim(fun.time.list)[1])
fun.time.list

fun.time.list$OTU_id <- ifelse(is.na(fun.time.list$Genus),ifelse(is.na(fun.time.list$Family),paste0(fun.time.list$number,'_o_',fun.time.list$Order),paste0(fun.time.list$number,'_f_',fun.time.list$Family)),paste0(fun.time.list$number,'_',fun.time.list$Genus))
fun.time.list$OTU_id

bac.time.list
fun.time.list

otu.list.time <- rbind(bac.time.list, fun.time.list)
dim(otu.list.time)
# write.xlsx(otu.list,'otu_id_book.xlsx')


OTU_id.list.time <- rbind(bac.time.list[c('OTU','OTU_id')],fun.time.list[c('OTU','OTU_id')])
OTU_id.list.time$OTU_id


## individual seed
bac.clean.ss
fun.clean.ss

summarize_phyloseq(bac.clean.ss)
bacASV<-otu_table(bac.clean.ss)
b.read.tab<-melt(colSums(bacASV))
max(colSums(bacASV))
min(colSums(bacASV))


bacASV[bacASV > 0]<- 1
colSums(bacASV)
b.ASVnum.tab<-melt(colSums(bacASV))
max(colSums(bacASV))
min(colSums(bacASV))

names(b.read.tab)[1] <- "Read Count"
names(b.ASVnum.tab)[1] <- "Number of ASVs"

b.summ.tab<-merge(b.ASVnum.tab, b.read.tab, by = "row.names")
names(b.summ.tab)[1] <- "Seed Sample"

b.summ.tab$`Seed Sample` <- factor(b.summ.tab$`Seed Sample`, levels = order.sample)

write.csv(b.summ.tab, "220520_Supplementary Table 3_bacteria.csv")




summarize_phyloseq(fun.clean.ss)
funASV<-otu_table(fun.clean.ss)
f.read.tab<-melt(colSums(funASV))
max(colSums(funASV))
min(colSums(funASV))


funASV[funASV > 0]<- 1
colSums(funASV)
f.ASVnum.tab<-melt(colSums(funASV))
max(colSums(funASV))
min(colSums(funASV))

names(f.read.tab)[1] <- "Read Count"
names(f.ASVnum.tab)[1] <- "Number of ASVs"

f.summ.tab<-merge(f.ASVnum.tab, f.read.tab, by = "row.names")
names(f.summ.tab)[1] <- "Seed Sample"

f.summ.tab$`Seed Sample` <- factor(f.summ.tab$`Seed Sample`, levels = order.sample)

write.csv(f.summ.tab, "220520_Supplementary Table 3_fungi.csv")


### pooled seed
summarize_phyloseq(bac.time.clean.ss)
bacASV<-otu_table(bac.time.clean.ss)
b.read.tab<-melt(colSums(bacASV))
max(colSums(bacASV))
min(colSums(bacASV))


bacASV[bacASV > 0]<- 1
colSums(bacASV)
b.ASVnum.tab<-melt(colSums(bacASV))
max(colSums(bacASV))
min(colSums(bacASV))

names(b.read.tab)[1] <- "Read Count"
names(b.ASVnum.tab)[1] <- "Number of ASVs"

b.summ.tab<-merge(b.ASVnum.tab, b.read.tab, by = "row.names")
names(b.summ.tab)[1] <- "Seed Sample"

write.csv(b.summ.tab, "220520_Supplementary Table 3_bacteria_pooled.csv")




summarize_phyloseq(fun.time.clean.ss)
funASV<-otu_table(fun.time.clean.ss)
f.read.tab<-melt(colSums(funASV))
max(colSums(funASV))
min(colSums(funASV))


funASV[funASV > 0]<- 1
colSums(funASV)
f.ASVnum.tab<-melt(colSums(funASV))
max(colSums(funASV))
min(colSums(funASV))

names(f.read.tab)[1] <- "Read Count"
names(f.ASVnum.tab)[1] <- "Number of ASVs"

f.summ.tab<-merge(f.ASVnum.tab, f.read.tab, by = "row.names")
names(f.summ.tab)[1] <- "Seed Sample"

write.csv(f.summ.tab, "220520_Supplementary Table 3_fungi_pooled.csv")
