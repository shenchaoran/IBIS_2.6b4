let Promise = require('bluebird')
let fs = Promise.promisifyAll(require('fs'))
let path = require('path')

let srcPath = path.join(__dirname, '..', 'src')
fs.readdirAsync(srcPath)
    .then(files => {
        return Promise.map(files, file => {
            if(file.indexOf('.h') !== -1) {
                let moduleName = file.replace(/.h/, '')
                let fname = file.replace(/.h/, '.module.f')
                let str = `      module ${moduleName}
        implicit none
        include '${file}'
        
      end module ${moduleName}`
                return fs.writeFileAsync(path.join(__dirname, '..', 'src', fname), str)
                    .then(() => {
                        console.log(`witefile succeed: ${fname}`)
                        return Promise.resolve()
                    })
                    .catch(e => {
                        console.log(`writefile failed: ${fname}`)
                        return Promise.resolve()
                    })
            }
            else 
                return Promise.resolve()
        }, {
            concurrency: 20
        })
        .then(() => {
            console.log('finished')
        })
        .catch(console.log)
    })
    .then(() => {
        console.log('finished')
    })
    .catch(console.log)