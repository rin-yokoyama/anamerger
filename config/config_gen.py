import argparse
import yaml

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-i','--input_file')
    parser.add_argument('-o','--output_file')
    parser.add_argument('-b','--block_to_replace')
    parser.add_argument('-k','--key_to_replace')
    parser.add_argument('-v','--value_to_replace')
    
    args = parser.parse_args()
    
    input_file = open(args.input_file)
    yaml_data = yaml.safe_load(input_file)
    input_file.close()
    
    yaml_data[args.block_to_replace][args.key_to_replace] = args.value_to_replace
    
    output_file = open(args.output_file,'w')
    output_file.write(yaml.dump(yaml_data,default_flow_style=False))
    output_file.close()

