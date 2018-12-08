import argparse
import pandas as pd

def parse_args():
    parser=argparse.ArgumentParser(description="aggregate GO terms across different Culture comparisons")
    parser.add_argument("--folders",nargs="+")
    parser.add_argument("--out_prefix")
    return parser.parse_args()

def main():
    args=parse_args()
    bp_dict=dict()
    mf_dict=dict()
    cc_dict=dict()

    bp_terms=set([])
    mf_terms=set([])
    cc_terms=set([])
    
    for dirname in args.folders:
        #GO biological process
        bp_data=pd.read_table('/'.join([dirname,'biological_process.txt']),header=0,sep='\t')
        for index,row in bp_data.iterrows():
            term=row['Term']
            bp_terms.add(term)
            logP=-1*row['logP']
            if term not in bp_dict:
                bp_dict[term]=dict()
            bp_dict[term][dirname]=logP
        print(','.join([dirname,'BP']))
        
        #GO molecular function
        mf_data=pd.read_table('/'.join([dirname,'molecular_function.txt']),header=0,sep='\t')
        for index,row in mf_data.iterrows():
            term=row['Term']
            mf_terms.add(term)
            logP=-1*row['logP']
            if term not in mf_dict:
                mf_dict[term]=dict()
            mf_dict[term][dirname]=logP
        print(','.join([dirname,'MF']))
        
        
        #GO cellular component
        cc_data=pd.read_table('/'.join([dirname,'cellular_component.txt']),header=0,sep='\t') 
        for index,row in cc_data.iterrows():
            term=row['Term']
            cc_terms.add(term)
            logP=-1*row['logP']
            if term not in cc_dict:
                cc_dict[term]=dict()
            cc_dict[term][dirname]=logP
        print(','.join([dirname,'CC']))

    #Write outputs
    bp_output=open(args.out_prefix+'.bp.txt','w')
    mf_output=open(args.out_prefix+'.mf.txt','w')
    cc_output=open(args.out_prefix+'.cc.txt','w')
    outputs=[bp_output,mf_output,cc_output]
    dirs=args.folders
    for output in outputs:
        output.write('Term\t')
        for dirname in dirs:
            output.write('\t'+dirname)
        output.write('\n')
    for term in bp_terms:
        bp_output.write(str(term))
        for dirname in dirs:
            if dirname in bp_dict[term]:
                bp_output.write('\t'+str(bp_dict[term][dirname]))
            else:
                bp_output.write('\tNA')
        bp_output.write('\n')
    for term in mf_terms:
        mf_output.write(str(term))
        for dirname in dirs:
            if dirname in mf_dict[term]:
                mf_output.write('\t'+str(mf_dict[term][dirname]))
            else:
                mf_output.write('\tNA')
        mf_output.write('\n')
    for term in cc_terms:
        cc_output.write(str(term))
        for dirname in dirs:
            if dirname in cc_dict[term]:
                cc_output.write('\t'+str(cc_dict[term][dirname]))
            else:
                cc_output.write('\tNA')
        cc_output.write('\n')
    
        

if __name__=="__main__":
    main()
    
