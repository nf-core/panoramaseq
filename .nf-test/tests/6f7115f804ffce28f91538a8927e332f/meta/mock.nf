import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies


// include test process
include { tsv_to_h5ad } from '/kyukon/scratch/gent/vo/000/gvo00027/vsc44685/panoramaseq/nf-core-panoramaseq/modules/local/anndata/make_h5ad/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    

    // process mapping
    def input = []
    
                input[0] = [
                    [ id: 'multi_h5ad' ],
                    [
                        file(params.panoramaseq_testdata_base_path + 'counts/S09_counts.tsv.gz', checkIfExists: true),
                        file(params.panoramaseq_testdata_base_path + 'counts/S10_counts.tsv.gz', checkIfExists: true)
                    ],
                    file(params.panoramaseq_testdata_base_path + 'barcodes_coords.csv', checkIfExists: true)
                ]
                
    //----

    //run process
    tsv_to_h5ad(*input)

    if (tsv_to_h5ad.output){

        // consumes all named output channels and stores items in a json file
        for (def name in tsv_to_h5ad.out.getNames()) {
            serializeChannel(name, tsv_to_h5ad.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = tsv_to_h5ad.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
  
}

def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
